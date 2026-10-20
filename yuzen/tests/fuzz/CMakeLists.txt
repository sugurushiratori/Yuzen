# ---------------------------------------------------------------------------
#  libFuzzer harnesses of the inputs (YUZEN_BUILD_FUZZERS=ON, clang only)
#
#    fuzz_input_toml  the .toml input: toml++ loader, Input::Read, CheckParameters
#    fuzz_input_inp   the same from a block-structured .inp (the ConfigBlock lexer/parser)
#    fuzz_restart     the restart readers: field*.dat, .vts, field*.state
#
#  Configure a build of its own, e.g.
#    CC=clang-18 CXX=clang++-18 cmake -S . -B build-fuzz -G Ninja \
#      -DCMAKE_BUILD_TYPE=RelWithDebInfo -DYUZEN_BUILD_FUZZERS=ON \
#      -DYUZEN_ENABLE_OPENMP=OFF -DYUZEN_WITH_AMGCL=OFF
#  write the seed corpora and the dictionaries with
#    python3 tools/fuzz_corpus.py --out fz --exe build/bin/yuzen
#  and run, e.g.
#    YUZEN_FUZZ_TMP=/dev/shm ASAN_OPTIONS=detect_leaks=0 \
#      build-fuzz/bin/fuzz_input_toml -detect_leaks=0 -close_fd_mask=3 \
#      -dict=fz/input.dict CORPUS_DIR fz/toml
#  (fuzz_restart: -dict=fz/restart.dict, seeds fz/restart). detect_leaks=0
#  in both places: a refused input leaves the frames it exited from
#  unreleased (see fuzz_common.h), and without the ASAN_OPTIONS one the leak
#  check at exit reports them, so that a run with no finding ends with exit
#  code 1 and an empty crash file.
#
#  Copyright 2017-2026 Suguru Shiratori, Tokyo City University
#  SPDX-License-Identifier: Apache-2.0
# ---------------------------------------------------------------------------

foreach(t fuzz_input_toml fuzz_input_inp fuzz_restart)
  if(t STREQUAL "fuzz_restart")
    add_executable(${t} fuzz_restart.cpp fuzz_exit.cpp)
  else()
    add_executable(${t} fuzz_input.cpp fuzz_exit.cpp)
  endif()
  target_link_libraries(${t} PRIVATE yuzencore)
  target_link_options(${t} PRIVATE -fsanitize=fuzzer -Wl,--wrap=exit)
endforeach()
target_compile_definitions(fuzz_input_toml PRIVATE YZFUZZ_EXT=".toml")
target_compile_definitions(fuzz_input_inp  PRIVATE YZFUZZ_EXT=".inp")
