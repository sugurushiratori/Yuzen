# Changelog

Every release of Yuzen, newest first. The manual's *Revision history*
appendix is the authoritative record and carries the same entries; the rule
for writing one (one line per user-visible change, and for any change that
moved a number, the case and the size of the difference) is set out there.

The section of a release is used as the text of its GitHub release, so it is
written for someone deciding whether to download it.

## [1.0.0]

The first public release.

Yuzen solves the long-wave (lubrication) equation for the height of a thin
liquid film on a substrate, coupled to a transport equation for the solute
concentration, in one and two horizontal dimensions.

### What it does

- **Models.** Well-mixed concentration in 1-D and 2-D (`Physics.ConcField.Model
  = "WM"`), and `"CZ"`, which resolves the concentration through the film
  thickness, in 1-D and 2-D.
- **Physics.** Capillary pressure, disjoining pressure (contact line and
  precursor film, calibrated to the equilibrium contact angle given in the
  input), vertical and in-plane gravity, centrifugal force (spin coating),
  soluto-capillary (Marangoni) stress, evaporation with a
  concentration-dependent viscosity, Navier slip at the substrate
  (`Physics.Slip`, well-mixed model), substrate topography, and a moving
  nozzle source that either delivers liquid at the concentration it lands on
  or at a concentration you give (`Physics.Source.ConcType`). Each term is
  switched on by its own `[Physics.*]` table; a table left out is a term
  switched off.
- **Boundaries.** `Symmetry`, `Periodic`, `Pinning`, `FreshInlet` and
  `Outlet`.
- **Time integration.** Crank–Nicolson, with the film-height equation driven to
  the implicit solution by a Newton–Kantorovich iteration with an analytic
  Jacobian for every term. The concentration is solved inside the same loop,
  so the converged state is the fully implicit coupled solution. The well-mixed
  solute balance is conservative to rounding. The capillary pressure, the
  disjoining pressure and the vertical gravity share one positive face
  mobility (Diez and Kondic), so the flux of a face is that mobility times the
  difference of the driving potential: the discrete free energy never rises,
  and a steady film is the minimum of it. The time step is optionally
  adaptive, from a local truncation-error estimate; a step whose iteration
  stalls far from convergence is cut and retried.
- **Convection schemes**, chosen per equation in `[Numerics.Convection]`:
  `UTOPIA` (default for the concentration), `Koren` (default for the film
  height), `KK`, `VanLeer`, `QUICK`, `Central4`, `Upwind1`, `Upwind2`,
  `Central2`.
- **Linear solvers.** Eigen (BiCGStab with block-Jacobi ILUT, always built)
  and AMGCL (optional, smoothed-aggregation AMG), chosen per equation. The
  concentration tries a light preconditioner first and falls back; the
  tolerance of each Newton pass follows the iteration (Eisenstat-Walker), and
  Newton starts from the film extrapolated from the last two steps.
- **Parallel.** Shared memory (OpenMP) only. With the Eigen backend the results
  are bit-identical for any thread count, and `ctest` asserts it by re-running
  the suite at 1 and 4 threads. A GCC build links LLVM's `libomp` when it is
  installed, and linear solves then run on min(T, n/500) threads
  (`YUZEN_GRAIN`): up to 56 threads no measured case got slower by more than
  7 %. A run uses at most the cores of one socket (with GCC's `libgomp`, also
  at most max(8, cells/1000)); asked for more with `OMP_NUM_THREADS`, it stops
  and says what to ask for (`YUZEN_THREAD_CAP=0` lifts the limit).
  `YUZEN_PROFILE=1` reports where the time of a run went, and `YUZEN_OPT`
  turns the algorithms above off for comparison.
- **Input and output.** TOML input; a misspelled or misplaced key stops the run
  with its full path and the likely intent, and `yuzen --check-input` checks a
  file without running it. Output in VTK (ParaView) and tab-separated text,
  with `History.dat`, `CalcInfo.dat`, probes and `JobInfo.txt`. The field is
  written every `TimeInterval` or at listed times (`Times`, or a `LogSchedule`
  spaced logarithmically), which cap only the step that would pass them, and
  with `ExactTimes` are landed on; `WritePotential` adds the driving
  potential (column `PHI`), and `Output.Energy` writes the free energy, its
  dissipation and their discrete balance (`Energy.dat`). Any field
  output is also a restart point, and a run stopped by `SIGINT`/`SIGTERM` and
  continued is bit-identical to one that was not stopped. Files are written
  under a temporary name and renamed when complete.

### Requirements

A C++17 compiler and CMake ≥ 3.16 (3.20 for `ctest --test-dir`). Eigen is
fetched at configure time if it is not installed; AMGCL is optional. The test
entries that run the solver need Python 3.

### Documentation

The manual is attached to this release in two editions built from the same
sources: `yuzen-manual.pdf` (A4) and `yuzen-manual.html` (one self-contained
file with search), each in English and in Japanese (`yuzen-manual-ja.pdf`,
`yuzen-manual-ja.html`; the English edition is the reference, and the build
checks that the two agree on every label, key, value, formula, listing and
number). The input-file reference covers every key, with its type, default
and range. `docs/theory.md`, `docs/numerics.md` and
`docs/validation.md` are shorter summaries.

### Verification

`ctest` runs 17 entries (about 25 minutes on two cores): unit tests that check
analytic solutions, conservation and every term's Jacobian against a finite
difference; the same problem posed twice in different ways (axes exchanged or
mirrored, 1-D against 2-D, other scales); problems with a closed-form answer;
manufactured solutions that converge at second order in space and time with
every term on; restart round trips; and input, error and signal handling.
The manual's *Verification* chapter gives every number with the command that
reproduces it, and its *What has not been verified* chapter lists the gaps.
The largest gaps are these:

- The contact line is verified in 1-D only.
- The order in time is measured with every term on, not for every combination
  of switches.
- The `amgcl` backend is exercised by no test, and it is not bit-reproducible
  across thread counts.
- `FreshInlet` and `Outlet` are checked against each other and against their
  balances, not against a known inflow.

### Citing

Please cite Yuzen using the metadata in `CITATION.cff`.
