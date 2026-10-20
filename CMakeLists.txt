#  Copyright 2017-2026 Suguru Shiratori, Tokyo City University
#  SPDX-License-Identifier: Apache-2.0

cmake_minimum_required(VERSION 3.16)
project(Yuzen LANGUAGES CXX)

# ---------------------------------------------------------------------------
#  yuzen : lubrication-approximation thin liquid film solver
#
#    input parser : src/IO/ConfigBlock (.inp) + toml++
#    utilities    : src/util/YuzenUtil (<filesystem>)
#    linear solver: Eigen / AMGCL
#    grid & fields: src/grid
# ---------------------------------------------------------------------------

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)
if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
  set(CMAKE_BUILD_TYPE Release CACHE STRING "" FORCE)
endif()
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)

#  MSVC reads source as the system code page unless told otherwise, and some
#  comments carry Greek letters and units (mu, superscripts). Without this it
#  warns on every one of them (C4819) and can mis-decode the file.
#  NOMINMAX keeps <windows.h> from defining min/max as macros.
if(MSVC)
  add_compile_options(/utf-8)
  add_compile_definitions(NOMINMAX _CRT_SECURE_NO_WARNINGS)
endif()

#  No fused multiply-add contraction. Where the hardware has FMA (every
#  AArch64 target: Apple silicon, Graviton, Grace; x86 built with -march that
#  includes it), GCC and Clang otherwise fuse a*b + c on their own, and may
#  fuse the same source expression differently where it is inlined in two
#  places. A run restarted from its state files then no longer matches the
#  run made straight through, to the bit (tools/proptest.py, seed 2283, on
#  AArch64). With contraction off every operation is rounded as written, as
#  on x86-64 by default, so this changes no result there.
if(CMAKE_CXX_COMPILER_ID MATCHES "GNU|Clang")
  add_compile_options(-ffp-contract=off)
elseif(CMAKE_CXX_COMPILER_ID STREQUAL "IntelLLVM")
  add_compile_options(-fp-model=precise -ffp-contract=off)
endif()

option(YUZEN_ENABLE_OPENMP "Build with OpenMP"          ON)
option(YUZEN_WITH_AMGCL    "Build the AMGCL linear-solver backend" ON)
option(YUZEN_STRICT_WARNINGS "Enable the warning set that Makehead had disabled" ON)
option(YUZEN_BUILD_TESTS   "Build the numerics test suite" ON)
set(YUZEN_OPENMP_RUNTIME "" CACHE FILEPATH
    "OpenMP runtime to link instead of the compiler's own (LLVM libomp.so, Intel libiomp5.so); empty: see YUZEN_PREFER_LIBOMP")
option(YUZEN_PREFER_LIBOMP "GCC: link LLVM's libomp instead of libgomp when it is installed" ON)

# ---------------------------------------------------------------------------
#  Public libraries
# ---------------------------------------------------------------------------
# Eigen (header-only). Fetched if it is not already on the system.
#   Ubuntu/Debian : apt install libeigen3-dev
#   Windows       : vcpkg install eigen3   /  conan
find_package(Eigen3 3.3 QUIET NO_MODULE)
if(NOT TARGET Eigen3::Eigen)
  message(STATUS "Eigen3 not found; fetching it with FetchContent")
  include(FetchContent)
  #  Only the headers are used. SOURCE_SUBDIR names a directory with no
  #  CMakeLists.txt, so Eigen's own build (its tests, and install rules that
  #  would put Eigen into the prefix of `cmake --install`) is not added.
  FetchContent_Declare(eigen
    GIT_REPOSITORY https://gitlab.com/libeigen/eigen.git
    GIT_TAG        3.4.0
    GIT_SHALLOW    TRUE
    SOURCE_SUBDIR  yuzen-headers-only)
  FetchContent_MakeAvailable(eigen)
  add_library(Eigen3::Eigen INTERFACE IMPORTED)
  set_target_properties(Eigen3::Eigen PROPERTIES
    INTERFACE_INCLUDE_DIRECTORIES ${eigen_SOURCE_DIR})
endif()

# AMGCL (header-only). Its OpenMP backend parallelizes the preconditioner too.
# On two-dimensional cases it is two to three times faster than Eigen, and it
# scales with threads. AMGCL_NO_BOOST is set, so Boost is not required.
# If it cannot be fetched, we fall back to a build with Eigen alone.
if(YUZEN_WITH_AMGCL)
  find_path(AMGCL_INCLUDE_DIR amgcl/amg.hpp)
  if(NOT AMGCL_INCLUDE_DIR)
    #  Cloned with git directly: a failed FetchContent_Populate is a configure
    #  error, so without network the fallback below would never be reached.
    set(_amgcl_src ${CMAKE_BINARY_DIR}/_deps/amgcl-src)
    if(NOT EXISTS ${_amgcl_src}/amgcl/amg.hpp)
      find_package(Git QUIET)
      if(GIT_FOUND)
        message(STATUS "AMGCL not found; fetching it with git")
        file(REMOVE_RECURSE ${_amgcl_src})
        execute_process(
          COMMAND ${GIT_EXECUTABLE} clone --quiet --depth 1 --branch 1.4.4
                  https://github.com/ddemidov/amgcl.git ${_amgcl_src}
          RESULT_VARIABLE _amgcl_rc
          TIMEOUT 300)
      endif()
    endif()
    if(EXISTS ${_amgcl_src}/amgcl/amg.hpp)
      set(AMGCL_INCLUDE_DIR ${_amgcl_src})
    endif()
  endif()
  if(AMGCL_INCLUDE_DIR)
    message(STATUS "AMGCL: ${AMGCL_INCLUDE_DIR}")
  else()
    message(WARNING "Could not fetch AMGCL; building with the Eigen backend only")
    set(YUZEN_WITH_AMGCL OFF)
  endif()
endif()

#  The libFuzzer harnesses (yuzen/tests/fuzz) need the whole library built
#  with the fuzzer's coverage and the sanitizers, so this is a build of its
#  own; clang only.
option(YUZEN_BUILD_FUZZERS "Build the libFuzzer harnesses of the inputs (clang)" OFF)
if(YUZEN_BUILD_FUZZERS)
  if(NOT CMAKE_CXX_COMPILER_ID MATCHES "Clang")
    message(FATAL_ERROR "YUZEN_BUILD_FUZZERS needs clang (libFuzzer)")
  endif()
  add_compile_options(-fsanitize=fuzzer-no-link,address,undefined
                      -fno-sanitize-recover=undefined -fno-omit-frame-pointer)
  add_link_options(-fsanitize=address,undefined)
endif()

# ---------------------------------------------------------------------------
#  yuzen proper
# ---------------------------------------------------------------------------
set(S ${CMAKE_SOURCE_DIR}/yuzen/src)

set(YUZEN_CORE_SRC
  ${S}/Yuzen.cpp
  ${S}/Version.cpp
  ${S}/IO/CheckParameters.cpp
  ${S}/IO/ConfigBlock.cpp
  ${S}/util/YuzenUtil.cpp
  ${S}/IO/FieldInfo.cpp
  ${S}/IO/Input.cpp
  ${S}/IO/Output.cpp
  ${S}/IO/Probe.cpp
  ${S}/IO/Restart.cpp
  ${S}/IO/RestartState.cpp
  ${S}/model/BoundaryCondition.cpp
  ${S}/model/Evaporation.cpp
  ${S}/model/ManufacturedSource.cpp
  ${S}/model/FluidProperty.cpp
  ${S}/model/Grids.cpp
  ${S}/model/InitialCondition.cpp
  ${S}/model/LubFields.cpp
  ${S}/model/NonDimensionalNumbers.cpp
  ${S}/model/Substrate.cpp
  ${S}/model/SolverModel.cpp
  ${S}/model/well_mixed/SourceBase.cpp
  ${S}/model/well_mixed/WM1D/Source_1D.cpp
  ${S}/model/well_mixed/WM1D/HeightEquation_1D.cpp
  ${S}/model/well_mixed/WM1D/ConcEquation_1D.cpp
  ${S}/model/well_mixed/WM1D/LaplacePressure_1D.cpp
  ${S}/model/well_mixed/WM1D/WellMixed_Model_1D.cpp
  ${S}/model/well_mixed/WM2D/HeightEquation_2D.cpp
  ${S}/model/well_mixed/WM2D/ConcEquation_2D.cpp
  ${S}/model/well_mixed/WM2D/LaplacePressure_2D.cpp
  ${S}/model/well_mixed/WM2D/WellMixed_Model_2D.cpp
  ${S}/model/well_mixed/WM2D/Source.cpp
  ${S}/model/ConcZ/CZ2D/ConcZ_Model_2D.cpp
  ${S}/model/ConcZ/CZ2D/CZ_ConcEquation_2D.cpp
  ${S}/model/ConcZ/CZ2D/CZ_HeightEquation_2D.cpp
  ${S}/numerics/CRS_Matrix.cpp
  ${S}/numerics/FiniteDifferences.cpp
  ${S}/numerics/ImplicitSolver.cpp
  ${S}/numerics/LinearSolver.cpp
  ${S}/numerics/LinearSolverBackend.cpp
  ${S}/numerics/NewtonKantorovich.cpp
  ${S}/numerics/SemiImplicit.cpp
  ${S}/numerics/TimeSteps.cpp
  ${S}/commonFunctions/ShapeFunctions.cpp
  ${S}/commonFunctions/GenerateShape.cpp)

add_library(yuzencore STATIC ${YUZEN_CORE_SRC})
target_include_directories(yuzencore PUBLIC
  ${CMAKE_SOURCE_DIR}/third_party
  ${S} ${S}/IO ${S}/util ${S}/grid ${S}/model ${S}/commonFunctions
  ${S}/model/well_mixed ${S}/model/well_mixed/WM1D ${S}/model/well_mixed/WM2D
  ${S}/model/ConcZ/CZ2D ${S}/numerics ${S}/MyMath)
target_compile_definitions(yuzencore PUBLIC _DEF_YUZEN_)
target_link_libraries(yuzencore PUBLIC Eigen3::Eigen)

if(YUZEN_WITH_AMGCL)
  target_include_directories(yuzencore PRIVATE ${AMGCL_INCLUDE_DIR})
  target_compile_definitions(yuzencore PRIVATE YUZEN_WITH_AMGCL AMGCL_NO_BOOST)
endif()

#  These are GCC/Clang spellings; MSVC rejects them outright, so guard the
#  block by compiler rather than applying it everywhere.
if(YUZEN_STRICT_WARNINGS)
  if(MSVC)
    target_compile_options(yuzencore PRIVATE /W3)
  else()
    target_compile_options(yuzencore PRIVATE
      -Wall -Wextra
      -Wimplicit-fallthrough
      -Wno-unused-parameter
      -Wno-deprecated
      -Wno-long-long)
  endif()
endif()

#  OpenMP is wanted but not required. AppleClang ships without an OpenMP
#  runtime, and a bare clang++ on Linux needs libomp installed separately, so
#  demanding it turned a missing optional dependency into a configure error on
#  those platforms. Build serially instead and say so: the solver is correct
#  either way, since results are bit-identical whatever the thread count.
if(YUZEN_ENABLE_OPENMP)
  find_package(OpenMP)
  #  GCC's own runtime, libgomp, makes every parallel region cost more the
  #  more threads take part (about 8 us at 56 threads, against 2 us with
  #  LLVM's libomp or Intel's libiomp5 on a Xeon Max 9480 socket), and
  #  cannot resize a team cheaply (util/Grain.h). So a GCC build links
  #  libomp instead when one is installed (libomp provides the entry points
  #  GCC's code calls); YUZEN_PREFER_LIBOMP=OFF, or a runtime named in
  #  YUZEN_OPENMP_RUNTIME, decides otherwise. clang uses libomp anyway, and
  #  the Intel compiler libiomp5.
  if(OpenMP_CXX_FOUND AND NOT YUZEN_OPENMP_RUNTIME AND YUZEN_PREFER_LIBOMP
     AND CMAKE_CXX_COMPILER_ID STREQUAL "GNU" AND NOT WIN32 AND NOT APPLE)
    file(GLOB _llvm_libdirs /usr/lib/llvm-*/lib /usr/lib64/llvm*/lib /opt/homebrew/opt/libomp/lib)
    list(SORT _llvm_libdirs ORDER DESCENDING)
    find_library(YUZEN_LIBOMP NAMES omp libomp.so.5 HINTS ${_llvm_libdirs})
    if(YUZEN_LIBOMP)
      set(_omp_runtime "${YUZEN_LIBOMP}")
    else()
      message(STATUS "OpenMP: LLVM libomp not found; GCC's libgomp is used "
                     "(slower parallel regions; install libomp for a faster build)")
    endif()
  endif()
  if(YUZEN_OPENMP_RUNTIME)
    set(_omp_runtime "${YUZEN_OPENMP_RUNTIME}")
  endif()
  if(OpenMP_CXX_FOUND AND _omp_runtime)
    #  Compiled with the compiler's OpenMP, linked against another runtime
    #  that provides its entry points (LLVM libomp and Intel libiomp5 provide
    #  GCC's GOMP_* as well): the runtime's barriers and thread pool instead
    #  of the compiler's own
    if(NOT EXISTS "${_omp_runtime}")
      message(FATAL_ERROR "YUZEN_OPENMP_RUNTIME: ${_omp_runtime} does not exist")
    endif()
    separate_arguments(_omp_flags NATIVE_COMMAND "${OpenMP_CXX_FLAGS}")
    target_compile_options(yuzencore PUBLIC ${_omp_flags})
    get_filename_component(_omp_dir "${_omp_runtime}" DIRECTORY)
    target_link_libraries(yuzencore PUBLIC "${_omp_runtime}" "-Wl,-rpath,${_omp_dir}")
    message(STATUS "OpenMP: compiled with ${OpenMP_CXX_FLAGS}, linked against ${_omp_runtime}")
  elseif(OpenMP_CXX_FOUND)
    target_link_libraries(yuzencore PUBLIC OpenMP::OpenMP_CXX)
  else()
    message(STATUS "OpenMP not found; building without it (serial run)")
    set(YUZEN_ENABLE_OPENMP OFF)
  endif()
endif()

add_executable(yuzen ${S}/mainProg.cpp)
target_link_libraries(yuzen PRIVATE yuzencore)
#  Windows: UTF-8 as the code page of the process, so that a path with
#  non-ASCII characters (from the command line or the UTF-8 input file) names
#  the same file as on Linux; see yuzen.manifest. MSVC's linker embeds a
#  .manifest given as a source.
if(MSVC)
  target_sources(yuzen PRIVATE ${S}/yuzen.manifest)
endif()

#  cmake --install build [--prefix DIR]: the solver into bin/, and the
#  license and notices into share/doc/yuzen/. The tools/ scripts and the
#  cases are run from the source tree and are not installed.
include(GNUInstallDirs)
install(TARGETS yuzen RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR})
install(FILES LICENSE NOTICE THIRD_PARTY_NOTICES.md README.md
        DESTINATION ${CMAKE_INSTALL_DATAROOTDIR}/doc/yuzen)

# ---------------------------------------------------------------------------
#  Tests
# ---------------------------------------------------------------------------
if(YUZEN_BUILD_TESTS)
  enable_testing()
  add_subdirectory(yuzen/tests)
endif()
if(YUZEN_BUILD_FUZZERS)
  add_subdirectory(yuzen/tests/fuzz)
endif()
