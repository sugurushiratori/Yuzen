# ---------------------------------------------------------------------------
#  The test suite
#
#  Every entry asserts: on the error against analytic or manufactured
#  solutions, on conservation and other invariants, or on the agreement of
#  the same problem posed in two ways.
#
#  Copyright 2017-2026 Suguru Shiratori, Tokyo City University
#  SPDX-License-Identifier: Apache-2.0
# ---------------------------------------------------------------------------

add_library(doctest_main STATIC doctest_main.cpp)
target_include_directories(doctest_main PUBLIC ${CMAKE_SOURCE_DIR}/third_party)
# doctest 2.4.11 forward-declares std::tuple, which recent MSVC warns about
# (C5285). The warning is in the third-party header, not in our tests.
if(MSVC)
  target_compile_options(doctest_main PUBLIC /wd5285)
endif()

set(YUZEN_TEST_SRC
  test_vector.cpp
  test_crs_matrix.cpp
  test_linear_solver.cpp
  test_timesteps.cpp
  test_finite_differences.cpp
  test_config.cpp
  test_boundary_condition.cpp
  test_grid.cpp
  test_newton_convergence.cpp
  test_mixed_velocity.cpp
  test_disjoining.cpp
  test_positivity_mobility.cpp
  test_restart.cpp
  test_source.cpp
  test_jacobian.cpp
  test_thread_cap.cpp
  test_opt.cpp)

add_executable(yuzen_tests ${YUZEN_TEST_SRC})
target_link_libraries(yuzen_tests PRIVATE doctest_main yuzencore)
target_include_directories(yuzen_tests PRIVATE ${CMAKE_SOURCE_DIR}/third_party)

add_test(NAME yuzen_tests COMMAND yuzen_tests)

# ---------------------------------------------------------------------------
#  Regression tests that actually run the solver
#
#  In the limit of a concentration that is uniform across the film
#  thickness, ConcZ solves the same problem as well-mixed. Both are in
#  conservative form, but the mobility and the flux are built independently
#  of each other, so the fact that they agree supports the validity of both.
#  If this stops matching, one of the two discretizations has been broken.
#
#  The tolerance 1e-5 is about 10 times the measured difference
#  (5.7e-07 .. 9.4e-07; it depends on Nx but not on Nz or on the thread
#  count). That difference comes from the significant digits of the output,
#  so anything far above it means the physics or the discretization has
#  changed.
#  Volume is conserved exactly by both models, so we check that one at
#  1e-10.
# ---------------------------------------------------------------------------
find_package(Python3 COMPONENTS Interpreter QUIET)
if(Python3_Interpreter_FOUND)
  add_test(NAME concz_vs_wm
    COMMAND ${Python3_EXECUTABLE}
            ${CMAKE_SOURCE_DIR}/tools/concz_vs_wm.py
            --exe $<TARGET_FILE:yuzen>
            --workdir ${CMAKE_CURRENT_BINARY_DIR}/concz_vs_wm
            --tol 1e-5 --vol-tol 1e-10)
  #  This runs the solver twice, so the default 1500 seconds is sometimes
  #  not enough
  set_tests_properties(concz_vs_wm PROPERTIES TIMEOUT 900)

  #  The same comparison in pseudo 3-D (the film spreads over x-y).
  #  The x-y cross terms (thickness equation) and the y-zeta cross terms
  #  (concentration equation) come into play here for the first time, so
  #  the 1-D test is no substitute for this one.
  add_test(NAME concz_vs_wm_2d
    COMMAND ${Python3_EXECUTABLE}
            ${CMAKE_SOURCE_DIR}/tools/concz_vs_wm.py
            --exe $<TARGET_FILE:yuzen>
            --nx 24 --ny 24 --nz 11 --dt 2e-6
            --workdir ${CMAKE_CURRENT_BINARY_DIR}/concz_vs_wm_2d
            --tol 1e-5 --vol-tol 1e-10)
  set_tests_properties(concz_vs_wm_2d PROPERTIES TIMEOUT 900)

  #  A restart must start from the state it was given, at the time it was
  #  given, and number its files on from where the run it continues
  #  stopped -- for both models, in 1-D and 2-D, from the .vts and from the
  #  field*.dat. None of that is reached by the unit tests. The .vts each
  #  source run writes is also checked against its field*.dat: the header's
  #  offsets, the block lengths, and the state and velocities it carries.
  add_test(NAME restart_roundtrip
    COMMAND ${Python3_EXECUTABLE}
            ${CMAKE_SOURCE_DIR}/tools/restart_roundtrip.py
            --exe $<TARGET_FILE:yuzen>
            --workdir ${CMAKE_CURRENT_BINARY_DIR}/restart_roundtrip)
  set_tests_properties(restart_roundtrip PROPERTIES TIMEOUT 900)

  #  A run stopped and restarted again and again (from the .vts and the
  #  field*.dat in turn, into the same result directory) must be the run it
  #  would have been straight through, to the bit: field files, History,
  #  CalcInfo (dt, iterations), probes, the steady-state and dry-out stops and the step
  #  counts -- with the adaptive step, evaporation, the CZ corrector and
  #  coupling passes, limited schemes and the moving nozzle on. Without the
  #  field*.state beside the field file, within the digits of that file.
  add_test(NAME restart_chain
    COMMAND ${Python3_EXECUTABLE}
            ${CMAKE_SOURCE_DIR}/tools/restart_roundtrip.py --chain
            --exe $<TARGET_FILE:yuzen>
            --workdir ${CMAKE_CURRENT_BINARY_DIR}/restart_chain)
  set_tests_properties(restart_chain PROPERTIES TIMEOUT 600)

  #  Inputs that would otherwise run with a silently wrong answer -- a misspelled key,
  #  a key in the wrong table, a one-sided Periodic, a CZ edge the model
  #  cannot impose, a pinned wall facing a film of full thickness -- must be
  #  refused with the reason, and a correctly pinned film must still run;
  #  likewise a manufactured solution whose concentration leaves (0, 1), and
  #  a value outside a range the manual promises, a grid past the size limit,
  #  a .inp string without its closing quote. A file that leaves optional
  #  blocks out must run with their defaults, and a field read back with its
  #  concentration outside [0, 1] must stop the run as diverged.
  add_test(NAME input_guards
    COMMAND ${Python3_EXECUTABLE}
            ${CMAKE_SOURCE_DIR}/tools/input_guards.py
            --exe $<TARGET_FILE:yuzen>
            --workdir ${CMAKE_CURRENT_BINARY_DIR}/input_guards)
  set_tests_properties(input_guards PROPERTIES TIMEOUT 600)

  #  What a new user, a full disk or a job scheduler does to a run: a
  #  command line, an input file or a result directory that is wrong, a
  #  restart file missing or cut short, a file size limit or a full device
  #  while writing, SIGINT / SIGTERM / SIGKILL in the middle of a run (the
  #  continued run must be the straight one to the bit), paths with spaces
  #  and Japanese characters. Each must stop with the path and the reason
  #  and leave no file cut short. About 10 s.
  add_test(NAME error_paths
    COMMAND ${Python3_EXECUTABLE}
            ${CMAKE_SOURCE_DIR}/tools/error_paths.py
            --exe $<TARGET_FILE:yuzen>
            --workdir ${CMAKE_CURRENT_BINARY_DIR}/error_paths)
  set_tests_properties(error_paths PROPERTIES TIMEOUT 600)

  #  The documentation of the input keys against what the parser does:
  #  every block and key of Input.cpp -- type, default, required or not,
  #  accepted values, range -- against the key tables of the manual, its
  #  counts, every \key{} it names, its TOML listings and the shipped cases,
  #  which the binary must also accept with --check-input. Parsing only; it
  #  takes about a second.
  add_test(NAME doc_keys
    COMMAND ${Python3_EXECUTABLE}
            ${CMAKE_SOURCE_DIR}/tools/doc_keys.py
            --exe $<TARGET_FILE:yuzen>
            --workdir ${CMAKE_CURRENT_BINARY_DIR}/doc_keys)
  set_tests_properties(doc_keys PROPERTIES TIMEOUT 300)
endif()

# The same problem posed twice must give the same answer: x and y
# exchanged, x mirrored, 1-D against a 2-D run uniform in y, a periodic film
# moved round the domain, other scales for the dimensionless form -- for
# both models -- two known answers (the body-force velocities on a flat
# film, the solute of a drying ConcZ film), ConcZ's advection against
# well-mixed's, and a retried time step against the same step taken directly.
if(Python3_FOUND)
  add_test(NAME equivalence
    COMMAND ${Python3_EXECUTABLE}
            ${CMAKE_SOURCE_DIR}/tools/equivalence.py
            --exe $<TARGET_FILE:yuzen>
            --workdir ${CMAKE_CURRENT_BINARY_DIR}/equivalence)
  set_tests_properties(equivalence PROPERTIES TIMEOUT 900)
endif()

# Problems with a closed-form answer: capillary leveling and its order in
# dx and dt, gravity, the traveling wave of in-plane gravity, spin-off,
# drying (well-mixed and ConcZ), and ConcZ diffusion through the thickness.
if(Python3_FOUND)
  add_test(NAME analytic
    COMMAND ${Python3_EXECUTABLE}
            ${CMAKE_SOURCE_DIR}/tools/analytic.py
            --exe $<TARGET_FILE:yuzen>
            --workdir ${CMAKE_CURRENT_BINARY_DIR}/analytic)
  set_tests_properties(analytic PROPERTIES TIMEOUT 900)
endif()

# The method of manufactured solutions, at an amplitude of 30 % of the film:
# second order in dx and in dt, for both models in 1-D and 2-D -- with the
# capillary pressure and both gravities, and then with every term of each
# model at once (and, in the well-mixed model, a concentration that varies),
# and ConcZ with a concentration that varies in the plane, mu(c), sigma(c),
# D(c) and Marangoni, its velocity profile checked as well.
if(Python3_FOUND)
  add_test(NAME mms
    COMMAND ${Python3_EXECUTABLE}
            ${CMAKE_SOURCE_DIR}/tools/mms.py
            --exe $<TARGET_FILE:yuzen>
            --workdir ${CMAKE_CURRENT_BINARY_DIR}/mms)
  set_tests_properties(mms PROPERTIES TIMEOUT 3600)
endif()

# The advection of the concentration by each convection scheme
# (Numerics.Convection): a front carried round a periodic film must stay
# within its bounds with the limited schemes, and a sine must converge at
# each scheme's order; an oblique sine carried and diffused on a 2-D film
# must converge at second order, settled, with UTOPIA and with Koren.
if(Python3_FOUND)
  add_test(NAME advection
    COMMAND ${Python3_EXECUTABLE}
            ${CMAKE_SOURCE_DIR}/tools/advection.py
            --exe $<TARGET_FILE:yuzen>
            --workdir ${CMAKE_CURRENT_BINARY_DIR}/advection)
  set_tests_properties(advection PROPERTIES TIMEOUT 600)
endif()

# The features no other test switches on: probes, the steady-state, dry-out
# and MaxSteps stops, and the growth of an adaptive time step.
if(Python3_FOUND)
  add_test(NAME features
    COMMAND ${Python3_EXECUTABLE}
            ${CMAKE_SOURCE_DIR}/tools/features.py
            --exe $<TARGET_FILE:yuzen>
            --workdir ${CMAKE_CURRENT_BINARY_DIR}/features)
  set_tests_properties(features PROPERTIES TIMEOUT 900)

  # The manufactured-solution source is generated with SymPy
  # (tools/mms_generate.py). Where SymPy is installed, check that the
  # generated file is what the generator makes of the equations now.
  execute_process(COMMAND ${Python3_EXECUTABLE} -c "import sympy"
                  RESULT_VARIABLE _yz_no_sympy OUTPUT_QUIET ERROR_QUIET)
  if(_yz_no_sympy EQUAL 0)
    add_test(NAME mms_source_current
      COMMAND ${Python3_EXECUTABLE} ${CMAKE_SOURCE_DIR}/tools/mms_generate.py --check)
  endif()
endif()

# Random valid inputs (tools/proptest.py): a fixed slice of 20 seeds -- both
# models, 1-D and 2-D, random terms, edges, shapes, schemes and steps --
# each run held to its properties: it ends or stops with a reason, nothing
# is NaN, the film stays positive, a bounded scheme keeps c in [0, 1], the
# volume and the solute are conserved, through open edges too; for some seeds a
# mirrored or transposed run agrees, two threads give the same bits, and
# the run restarted twice is the run straight through.
if(Python3_FOUND)
  add_test(NAME proptest
    COMMAND ${Python3_EXECUTABLE}
            ${CMAKE_SOURCE_DIR}/tools/proptest.py --quick
            --exe $<TARGET_FILE:yuzen>
            --workdir ${CMAKE_CURRENT_BINARY_DIR}/proptest)
  set_tests_properties(proptest PROPERTIES TIMEOUT 600)
endif()

# The solver runs above make thousands of short OpenMP loops. The default
# wait policy spins at every barrier, and once the machine has fewer free
# cores than threads -- ctest -j, a second job, a busy CI runner -- each
# barrier waits a scheduler slice instead: concz_vs_wm takes 246 s instead of
# 5 s with two threads on a machine with one core busy, and equivalence and
# analytic run into their timeouts. Passive waiting costs little when the
# cores are free, and the answers do not depend on it.
if(Python3_FOUND)
  foreach(t concz_vs_wm concz_vs_wm_2d restart_roundtrip restart_chain input_guards
            equivalence analytic mms features advection proptest)
    if(TEST ${t})
      set_property(TEST ${t} APPEND PROPERTY ENVIRONMENT "OMP_WAIT_POLICY=PASSIVE")
    endif()
  endforeach()
endif()

# Confirm in separate processes that changing the thread count does not
# change the results
if(YUZEN_ENABLE_OPENMP)
  add_test(NAME yuzen_tests_1thread COMMAND yuzen_tests)
  set_tests_properties(yuzen_tests_1thread PROPERTIES ENVIRONMENT "OMP_NUM_THREADS=1")
  add_test(NAME yuzen_tests_4thread COMMAND yuzen_tests)
  set_tests_properties(yuzen_tests_4thread PROPERTIES ENVIRONMENT "OMP_NUM_THREADS=4")
endif()
