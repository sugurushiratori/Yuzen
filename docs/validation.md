# Validation results

This document collects **what has been checked in Yuzen, how, and what the
numbers came out at**. It is a shorter companion to the manual's chapters
"Verification" and "What has not been verified" (`docs/manual/`), which give
every check in full; where the two differ, the manual is authoritative. The
governing equations and the scaling are in [`theory.md`](theory.md); the
discretization and the iterative methods are in [`numerics.md`](numerics.md).
Nothing is re-derived here: only the results and the procedures are given, and
every number can be reproduced from the tests, cases and tools in the
repository (§15).

---

## 1. How the claims are graded

Claims come in strengths, in descending order: agreement with an analytic
solution (§5.1, §9), order of convergence (§5.2, §9, §10), invariants and
identities (§4, §8), and regression fixing alone (input interpretation, output
format, time-step control), which records the current behavior and says
nothing about it being right.

**The strong claims.** That the disjoining-pressure coefficient `Dp` agrees
with the equilibrium contact angle asked for has been confirmed
**analytically** over 36 parameter combinations (§5.1). For the discretization
of capillary plus disjoining pressure, the residual velocity on an analytic
equilibrium droplet **falls to 0 at order 1.94** (§5.2) — if the coefficient,
the sign or the discretization were wrong in any way, that residual would not
vanish under mesh refinement. The reductions and BiCGStab return
**bit-identical values under a change of thread count** (§4.2), tested with
`==`.

**The weak claims.** For the Cox-Voinov law (§6) only the form of the relation
and the correlation are examined; the absolute value of the logarithmic factor
`ln(L/l)` has not been verified independently. The cost of the coupling (§7) is
a measurement of iteration counts over one 1-D, well-mixed case family. Input
interpretation, output format and time-step control are regression fixes and
have no external ground truth. What has not been verified is listed in §16.

---

## 2. The test suite

```sh
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build -j
ctest --test-dir build
```

`ctest` registers **seventeen entries** (sixteen where SymPy is not installed;
all but the three test-binary entries need Python).

| Entry | What it requires |
|---|---|
| `yuzen_tests`, `yuzen_tests_1thread`, `yuzen_tests_4thread` | the unit tests below, in three separate processes — the environment's thread count, `OMP_NUM_THREADS=1` and `=4` — because **invariance under a change of thread count is itself one of the claims** (§4.2) |
| `concz_vs_wm`, `concz_vs_wm_2d` | the two concentration models agree where the concentration cannot affect the film (§8) |
| `restart_roundtrip` | a restart of either model, 1-D and 2-D, from the `.vts` (in `CZ` also one without the thickness and substrate) and from `field*.dat` starts from the state and time it was given, and each `.vts`, read through its header offsets, holds what the field file of the same instant holds (§12) |
| `restart_chain` | a run stopped and restarted several times gives the same files, bit for bit, as one run straight through (§12) |
| `input_guards` | inputs that would otherwise run with a silently wrong answer are refused with the reason (§12) |
| `error_paths` | what a new user, a full disk or a job scheduler does to a run stops it with the path and the reason, leaves no file cut short, and an interrupted run continued is the straight run to the bit (§12) |
| `equivalence` | the same problem posed twice gives the same answer (§8) |
| `analytic` | problems with a closed-form answer get it, at the order the scheme promises (§9) |
| `mms` | manufactured solutions with every term on converge at second order in space and time (§10) |
| `advection` | the convection schemes keep a step bounded where they promise to and converge on a smooth profile at their order; the well-mixed concentration converges in 2-D at second order with convection and diffusion together (§11) |
| `features` | probes, the stop conditions, the adaptive step and its clock, edited restart files (§12) |
| `proptest` | random valid inputs end cleanly, keep the film positive and the volume and solute, and are invariant where they must be (§13) |
| `doc_keys` | runs no solver: every key the parser reads, with its type, default, range and accepted values, against the key tables of the manual and of [`numerics.md`](numerics.md), the manual's TOML listings and every file in `cases/`, each also read by `yuzen --check-input` |
| `mms_source_current` | only where SymPy is installed: the generated manufactured source is what the generator makes of the current equations |

Calling the test binary directly prints the counts.

```
$ ./build/bin/yuzen_tests
[doctest] test cases:    128 |    128 passed | 0 failed | 0 skipped
[doctest] assertions: 183297 | 183297 passed | 0 failed |
```

| File | Cases | Assertions | What it covers |
|---|---:|---:|---|
| `test_grid.cpp` | 9 | 6,252 | Grid indices, moves and coordinates; where fields sit |
| `test_linear_solver.cpp` | 6 | 1,663 | The analytic solution of the 1-D Poisson problem, 2nd-order convergence, bit reproducibility |
| `test_finite_differences.cpp` | 11 | 171,651 | Interpolation coefficients and gradients against analytic values; every convection scheme: exactness, bounds, symmetry, frozen weights |
| `test_boundary_condition.cpp` | 11 | 1,025 | Matrix consistency over every combination of boundary conditions; the ghost values they give, inlets on every edge, mirrors and transposes; the boundary rows satisfied exactly after a solve |
| `test_mixed_velocity.cpp` | 3 | 216 | That the mixed velocity is the plain sum of the individual terms |
| `test_disjoining.cpp` | 4 | 179 | Calibration of the disjoining pressure and the residual velocity on the equilibrium shape |
| `test_restart.cpp` | 12 | 193 | The correspondence between lines of the ASCII restart file and grid points; the marker search of the binary reader; the restart state file; damaged restart files refused |
| `test_config.cpp` | 11 | 76 | Interpretation of the input file, in TOML and in the block-structured `.inp` format, the report of unknown entries, and the strict conversions |
| `test_vector.cpp` | 5 | 34 | Basic operations and bit reproducibility of the reductions |
| `test_positivity_mobility.cpp` | 11 | 858 | That the face mobility is the positivity-preserving form of Diez & Kondic; with Navier slip, against a reference without cancellation, continuous at its series, its derivative; zero against a `Pinning` ghost |
| `test_crs_matrix.cpp` | 4 | 28 | Element access in CRS and `MatVec` |
| `test_timesteps.cpp` | 13 | 169 | Time advance and output intervals; the output schedule after a restart; the clock advances by the step solved |
| `test_newton_convergence.cpp` | 9 | 23 | How the outer iteration exits, and when a stagnation fails the step |
| `test_source.cpp` | 2 | 15 | The nozzle arm's angle against the integral of its rotation rate |
| `test_opt.cpp` | 2 | 29 | The `YUZEN_OPT` and `YUZEN_GRAIN` switches parsed |
| `test_thread_cap.cpp` | 5 | 65 | The limit on the thread count |
| `test_jacobian.cpp` | 10 | 821 | Each term's Jacobian against a finite difference of the term, 1-D and 2-D, with every convection scheme, with and without slip; no flux through a wall; the capillary and gravity `A h` face by face; each term's face fluxes against its `A h` |
| Total | 128 | 183,297 | |

The assertion count reflects the number of loop passes inside a case, not the
number of properties: `test_finite_differences.cpp` is 94 % of the assertions
because it checks every convection scheme at every face of many random
profiles, and nothing like 94 % of the coverage.

---

## 3. Grid and finite differences

`test_grid.cpp` builds a `RectGrid<3>` with `Nx, Ny, Nz = 9, 7, 5` and checks
the following at **every grid point of all 8 site types (Node / Xedge /
Yedge / Zedge / Xface / Yface / Zface / Cell)**.

* The linear index and the grid coordinates are inverses of one another.
* Neighbors close periodically (`Nx` moves of 1 cell return to the start).
* A 1-cell move and a half-cell move both return to the start when reversed;
  for a half-cell move the destination site type matches as well.
* The index moves under a half-cell move exactly as the rule says — "a move of
  +half a cell leaves the index alone when that axis's stagger bit is 0 and
  gives +1 when the bit is 1; a move of -half a cell does the reverse".
* The coordinates follow the staggered layout. Nothing degenerates in a
  direction that holds only one line, such as `Ny = 1`.
* Fields are initialized to 0 and can be read back by index. The components of
  a vector field sit in different places from one another.

`test_finite_differences.cpp` pins down the difference side. Over the 7 fixed
schemes (`Upwind_1st` / `Upwind_2nd` / `Central_2nd` / `Upwind_3rd` / `QUICK` /
`KK` / `Central_4th`) it checks that the interpolation coefficients sum to 1
(so that a constant is preserved), that the gradient of a linear function is
reproduced exactly, that the gradient of a quadratic agrees to 2nd order or
better, and that the `Node` variant is the average of the `Xedge` values on
either side of it. The two limited schemes, Koren and van Leer, are checked
with the others in §11.

---

## 4. Linear algebra and reproducibility

### 4.1 The linear solver (analytic solution and order of convergence)

`test_linear_solver.cpp` drives the `eigen` backend through the
`ILinearSolverBackend` interface on the 1-D Poisson problem `u'' = 2` on
`[-1, 1]` with `u(±1) = 0` (analytic solution `u(x) = x^2 - 1`). All 3
preconditioners (`None` / `Jacobi` / `ILU`) converge, with a
residual < 1e-11 and a maximum error against the analytic solution at the
interior points < 1e-4 (`n = 202`, `dx = 0.01`; 2.5e-5 measured). Refining
through `n = 52, 102, 202, 402`, the maximum error is **more than halved**
every time (2nd-order convergence). A zero right-hand side returns the trivial
solution; the test asserts nothing about the iteration count (`eigen` takes 520
iterations to it). `Solve()` restores the matrix and the right-hand side to
their state before the call (the diagonal scaling of the preconditioner
overwrites both).

### 4.2 Determinism and bit reproducibility

The properties pinned down by the unit tests, compared with `==`, not with
`Approx`.

| Subject | Condition | Claim |
|---|---|---|
| `Vector` inner product / `Sum` / `Norm` | `N = 100003` | **Exactly equal** for thread counts 1, 2, 4, 8 |
| Linear solver, `eigen` (ILU) | `n = 402` | For thread counts 1, 2, 4, 8, every component of the solution agrees, **and so do the iteration count and the residual** |
| Linear solver, `eigen` (ILU) | The same problem twice | Iteration count, residual and solution agree exactly |

The reductions fix the order of summation with a block decomposition that does
not depend on the thread count (at least 512 elements per block, at most 256
blocks). `N = 100003` is not a multiple of 256 and is chosen so as to straddle
block boundaries. Measured over the solver as a whole (`tools/compare.py`,
maximum difference over column scale):

| Comparison | Result |
|---|---|
| `leveling_1d_reg`, 1 thread vs 2 threads (`eigen`) | 0.000e+00 |
| `leveling_1d_reg`, run twice on 2 threads | 0.000e+00 |
| `allphys_2d_reg`, 1 thread vs 2 threads (`eigen`) | 0.000e+00 |

### 4.3 Backend agreement

`leveling_1d_reg` (Nx = 300, 20 output intervals) was solved with `eigen` and
with `amgcl` on one thread and the outputs compared.

| What is compared | `eigen` vs `amgcl` |
|---|---:|
| Physical quantities in `field*.dat` | 7.033e-11 (`Ux`) |
| `History.dat` | 9.615e-02 (`ResinVolErr`) |
| Time and step in `CalcInfo.dat` | 0.000e+00 |

The 7.033e-11 is the linear solver's tolerance (`Epsilon = 1e-13`), not an
`O(dt)` discrepancy. `VolErr`, `ResinVolErr` and the conservation-error columns
generally hover around zero, so their relative difference is magnified: nothing
flows out of this case, both are rounding (at most 7e-15) and the two backends
differ on them by 6e-16. `tools/compare.py` reports the diagnostic
columns, `BCG_iter` among them, but keeps them out of the pass/fail decision.

The same comparison ruled out an in-house BiCGStab backend tried during
development: it agreed with `eigen` as closely (1.220e-10 in `Ux`) but needed
**3676** BiCGStab iterations over the 20 output intervals against **21** for
`eigen`, and it was left out of the release.

**`amgcl` is not thread-invariant.** Its internal reductions are not
order-fixed, so changing the number of threads can change the low-order bits.
The same `leveling_1d_reg` on 1 thread and on 4 differs by **6.089e-11** in the
physical quantities (`Ux`) and **9.615e-2** in the error measures
(`ResinVolErr`, rounding here: at most 6e-15, the two runs differing by 6e-16).
On 2 threads it gave the same bits as on 1, and so did `allphys_2d_reg` on 1,
2 and 4: whether the bits move depends on how the library splits the work, and
nothing guarantees that they do not. For the identical comparison `eigen` gives
0.000e+00 at any thread count (§4.2). `ctest` passes on an AMGCL build as well, but these numbers
come from one case run by hand (§16).

### 4.4 Velocity composition

The mixed velocity is the film's velocity in the output, and the Courant check
is formed from it (the well-mixed concentration itself is carried by the face
volume fluxes, and what crosses an open edge is counted from them,
[`numerics.md`](numerics.md) §1.7), so if the rule by which the velocities written by the individual physical
terms are combined goes wrong, all of these are wrong. The sign
convention of the form `-div(...)` is already inside each term's `_Magnitude`,
so the composition must be a plain sum. `test_mixed_velocity.cpp` gives 6
terms distinct values (capillary pressure 1, Marangoni 2, vertical gravity 4,
in-plane gravity 8, centrifugal force 16, disjoining pressure 32) and requires
the **plain sum, 63**; vertical gravity alone set to `+3` gives `+3` and
centrifugal force alone set to `-5` gives `-5` (**signs pass through
unchanged**); the same holds in 2-D for both components.

### 4.5 Boundary conditions, restart files and input parsing

* `test_boundary_condition.cpp` — over every combination of the 5 boundary
  conditions (`Symmetry` / `Periodic` / `Pinning` / `FreshInlet` / `Outlet`)
  and the 3 fields (`Height` / `HeightJacobi` / `Conc`): the column indices
  stay within range, every ghost row carries exactly 1 equation, and a
  symmetric boundary condition gives a matrix that is its own mirror image. It
  also checks the values: a uniform film is returned at every ghost point,
  corners included; an inlet on each edge holds its face at the value of the
  field it was built from; mirroring or transposing the field with the edge
  types moves every ghost value with it; and the boundary rows are satisfied
  exactly after a solve. These found that an inlet on Xmax or Ymax was held at
  the mean of its two ghost layers rather than of the two points either side of
  its face, and that its face velocity was overwritten with the one a cell
  inside; both corrected (the difference from the mirrored run fell from 0.25 to
  2e-9).
* `test_restart.cpp` — `field*****.dat` is both an output file and the input
  format of an ASCII restart, so **the line number written and the grid point
  read back have to be the same point**: (a) in 1-D, line `i` lands at x index
  `i+1`; (b) an asymmetric shape written out and read back returns to the same
  places (1e-13 relative); (c) the decomposition of the line number is correct
  on a 2-D grid with **`Nx != Ny` (12 x 5)**, where a mix-up between dividing by
  `Nx` and by `Ny` shows. It also covers the marker search of the binary
  reader, the restart state file, and damaged restart files being refused.
* `test_config.cpp` — scalars, vectors, nested blocks and repeated blocks; a
  `.toml` input producing the same tree as the equivalent `.inp`; the report of unknown entries; and the
  **strict conversions**: a value is converted whole or refused, naming the key
  and the value (trailing text, `nan`/`inf`, fractional counts, a list for a
  scalar and dates are refused).

---

## 5. Disjoining pressure and the contact line

This is the core of the validation of the contact-line model. The tests are in
`yuzen/tests/test_disjoining.cpp`.

### 5.1 The calibration of the coefficient can be checked analytically

Equilibrium means uniform pressure, that is
$`H'' + \mathrm{Dp}\,\hat{\Pi}(H) = K`$
($`K`$ constant = the pressure, $`\hat{\Pi}(H) = (H_*/H)^{n} - (H_*/H)^{m}`$).
In the limit where the precursor film extends to infinity ($`K = 0`$),
multiplying by $`H'`$ integrates the equation once:

```math
\frac{1}{2} H'^{2} = -\mathrm{Dp}\int_{H_*}^{H} \hat{\Pi}\,ds
```

```math
\int_{H_*}^{\infty} \hat{\Pi}\,ds = \frac{H_*(m-n)}{(n-1)(m-1)} \; (< 0)
```

```math
\Rightarrow\quad (H')^{2} \to \frac{2\,\mathrm{Dp}\,H_*(n-m)}{(n-1)(m-1)}
  = \mathrm{As}\,\theta_e^{2}
```

Since the physical slope is the dimensionless slope $`H'`$ multiplied by
$`\varepsilon = h_0/L_0`$, **the slope on the outside of the film is exactly
$`\theta_e`$**. The coefficient the implementation uses,
$`\mathrm{Dp} = \mathrm{As}\,(n-1)(m-1) \,/\, (2 H_* (n-m)) \cdot \theta_e^{2}`$,
agrees with this identically.
The test sweeps the **36 combinations** $`\theta_e \in \{2, 5, 10, 20\}`$
degrees, $`H_* \in \{0.01, 0.05, 0.2\}`$ and
$`(n, m) \in \{(3,2), (4,3), (9,3)\}`$, and checks that $`\mathrm{Dp}`$ agrees
with the formula above to within 1e-14 relative and that
$`\varepsilon\sqrt{-2\,\mathrm{Dp}\int\hat{\Pi}}`$ agrees with $`\theta_e`$ in
radians to within **1e-12 relative**. That disabling the disjoining pressure
makes $`\mathrm{Dp}`$ exactly 0 is fixed by a separate case.

### 5.2 No flow arises on an equilibrium droplet (2nd-order convergence)

**How the droplet is built.** On the flat film $`H_\mathrm{f}`$ surrounding the
droplet, $`H' = H'' = 0`$, so $`K = \mathrm{Dp}\,\hat{\Pi}(H_\mathrm{f})`$. At
the crest $`H_\mathrm{max}`$, $`H' = 0`$ as well, so

```math
K (H_\mathrm{max} - H_\mathrm{f})
  - \mathrm{Dp}\left[ W(H_\mathrm{max}) - W(H_\mathrm{f}) \right] = 0,
\qquad W(H) = \int_{H_*}^{H} \hat{\Pi}\,ds
```

links $`H_\mathrm{f}`$ and $`H_\mathrm{max}`$. Give $`H_\mathrm{max}`$, fix
$`H_\mathrm{f}`$ by bisection, and integrate
$`H'' = K - \mathrm{Dp}\,\hat{\Pi}(H)`$ from the foot to the crest, and the
shape is determined. **Equilibrium droplets are a 1-parameter family** and can
be built from 1 ordinary differential equation.

**What is measured.** Put that shape on the mesh and add up the capillary and
disjoining-pressure velocities. In the continuous equilibrium they cancel, so
what is left is the truncation error of the discretization. The conditions are
`theta_e = 10` degrees, `H* = 0.02` (= 0.2 µm), `Hmax = 1.0` (= 10 µm), on a
domain of ±300 µm. Run the solver, read `Ux` out of `field00000.dat` and make it
dimensionless with the velocity scale `U0 = 5.6069e-08` m/s:

| Nx | dx [µm] | residual velocity (dimensionless) | ratio | order |
|---:|---:|---:|---:|---:|
| 1000 | 0.6 | 18.94 | | |
| 2000 | 0.3 | 5.074 | 3.73 | 1.90 |
| 4000 | 0.15 | 1.294 | 3.92 | 1.97 |

Taken from 1000 to 4000 as a whole, the **order is 1.94** (1.97 on the finer
pair). The residual
velocity is largest at the contact line (`H` ≈ 2`H*`), not at the crest and not
at the foot.

The same quantity is measured independently inside the test itself: the
droplet is constructed on the spot with **4th-order Runge-Kutta**, the
velocities of `LaplacePressure_WM1D` and `DisjoiningPressureKernel<WM1D_Grid>`
are added, and the maximum is taken (`YUZEN_DBG_EQ=1` prints the breakdown).

```
  [eq] nx=1000 umax=2.0295e+01  at X=-6.843746e-02 H=3.953799e-02
  [eq] nx=2000 umax=5.1789e+00  at X=6.869474e-02 H=3.599702e-02
  [eq] nx=4000 umax=1.3013e+00  at X=6.865186e-02 H=3.656075e-02
```

The ratios are 3.92 and 3.98. **The solver and the independent integration
inside the test agree to within 0.6 % at Nx = 4000, 2 % at Nx = 2000 and 7 %
at Nx = 1000**: the residual is the difference of a capillary and a
disjoining-pressure velocity each about 170 times as large, so the two
constructions' interpolation errors show in it. The test asserts that refining always reduces the
residual and that the ratio lies between 3.0 and 5.0 (the band is wide because
the interpolation error of the shape construction rides along); a wrong
coefficient or a wrong sign puts the ratio near 1.

### 5.3 It does not move under time evolution either

The same equilibrium shape (Nx = 2000) is advanced in time for 0.01 s.

| t [s] | max&#124;Δh&#124; [µm] | relative | residual velocity (dimensionless) |
|---:|---:|---:|---:|
| 0 | 0 | 0 | 5.074 |
| 2.5e-3 | 4.23e-5 | 4.23e-6 | 0.709 |
| 5.0e-3 | 5.36e-5 | 5.36e-6 | 0.418 |
| 7.5e-3 | 6.02e-5 | 6.02e-6 | 0.296 |
| 1.0e-2 | 6.48e-5 | 6.48e-6 | 0.229 |

The shape barely moves (6.5e-6 relative) and the residual velocity **decreases
monotonically**, as the solution settles from the continuous equilibrium
toward the **discrete** one, the film of uniform discrete potential (one face
mobility for the capillary and the disjoining pressure, [`numerics.md`](numerics.md)
§1.4). The apparent contact angle (§5.5) goes from `9.56613` degrees to
`9.56610` degrees — **unchanged to five digits**. With the
default steady-state test left enabled the run is declared steady almost
immediately, so it has to be switched off to produce this table (§15.3).

A uniform precursor film `H = H*` is an equilibrium in its own right
(`Pihat(H*) = 0`). The test fills the velocity field with 1.0 first, then calls
`CalVelocity` and checks that the disjoining pressure's velocity is **below
1e-12** at every interior `Xedge`.

### 5.4 Two traps in constructing the equilibrium droplet

Each one **destroys the convergence-order measurement itself**;
`test_disjoining.cpp` carries both as comments.

**(a) The shape must be interpolated with cubic Hermite.** The 2nd difference
divides by `dx^2`. The error of a linear interpolation is `O(dX^2)` in the step
`dX` of the ODE integration — **a fixed height error, independent of `dx`** —
so divided by `dx^2` it grows as `dx^-2` and **the residual velocity increases
as the mesh is refined**. Cubic Hermite makes it `O(dX^4)`.

**(b) The initial value at the foot must be tuned so that the crest lands
exactly on a grid point.** Scaling the foot's initial value `eta0` by `c`
shifts the shape along `X` by `-ln(c)/lambda`, which is the knob that puts the
last integration point on the crest (`H' = 0`). Skip it and the mirrored
profile has a kink at the crest, again an error independent of `dx`.

### 5.5 The apparent contact angle and the model's precursor-film bias

The equilibrium droplet of a model with a finite precursor film does not come
out at `theta_e` unless `Hmax/H*` is large enough. The slope of the outer
parabola extrapolated to `H = 0` — the contact angle measured in an experiment
— is, in closed form, `theta_cap = theta_e (1 - H*/(2 Hmax))`. The table
produced by `tools/make_equilibrium_1d.py --out /dev/null --sweep`
(`theta_e = 10` degrees, `H* = 0.02`, `n = 3`, `m = 2`):

| Hmax/H\* | 2 | 5 | 10 | 20 | 50 | 100 | 1000 |
|---|---:|---:|---:|---:|---:|---:|---:|
| apparent angle [deg] | 7.5000 | 9.0000 | 9.5000 | 9.7500 | 9.9000 | 9.9500 | 9.9950 |
| maximum slope of the shape [deg] | 0.6320 | 3.6012 | 5.4236 | 6.7634 | 7.9627 | 8.5652 | 9.5504 |

The numerically integrated apparent angle agrees with the closed form to 4
digits; the maximum slope converges far more slowly. **The rules of thumb when
setting up a case** ([`theory.md`](theory.md) §4.4), in terms of the transition
length $`\ell = h_*/\theta_e`$:

| For what | Rule of thumb |
|---|---|
| Resolving the contact line | $`\Delta x \le \ell/5`$ |
| Getting the angle | $`H_\mathrm{max}/h_* \ge 50`$ |

At $`\theta_e = 10`$ degrees and $`h_* = 0.2`$ µm, $`\ell = 1.15`$ µm. The case
of §5.2 (`h0 = 10` µm, `h* = 0.2` µm, domain 600 µm, Nx = 2000) sits exactly on
that boundary.

---

## 6. Spreading (Cox-Voinov)

### 6.1 What is measured, and the measurement

Tanner's law (`R ~ t^(1/7)` for a 2-D ridge) is a law of complete wetting
(`theta_e = 0`), and this form of disjoining pressure cannot represent it
because `theta_e -> 0` implies `Dp -> 0`. What can be measured in partial
wetting is the **Cox-Voinov relation** `theta_d^3 - theta_e^3 = 9 ln(L/l) * Ca`,
with `Ca = mu U / sigma`. Each output time of one run is one `(theta_d, U)`
point.

`cases/eq/spread_n500.toml` — an equilibrium droplet built for `theta_e = 20`
degrees is the initial condition and the run is made with `ContactAngle = 10.0`.
`tools/measure_spreading.py` fits a parabola by least squares to the outside of
the droplet at each time and extrapolates it to `h = 0`, giving the base radius
`R` and the apparent angle `theta`.

```
  theta : 19.1201 -> 9.4397 deg     (t = 0 .. 1 s)
  R     : 59.744  -> 85.507 um
  volume proxy theta*R^2 : 1.1911e-09 -> 1.2046e-09   (conserved to within 1.1 %)

  theta^3 - theta_e^3 = 19.6272 * Ca + 5.97e-05
    => 9 ln(L/l) = 19.627  =>  ln(L/l) = 2.181   (L/l = 8.9)
  correlation = 0.99791   (108 points, Ca = 1.41e-05 .. 1.12e-03)
```

**A correlation of 0.998 over two decades of Ca.** The intercept 5.97e-05 is
about 1 % of the scale of `theta^3` (9.44 degrees = 0.165 rad, giving 4.5e-3),
which shows that the equilibrium angle has been taken correctly. The factor
`ln(L/l) = 2.18` is of the same order as the naive estimate
`ln(85.51 µm / 1.196 µm) = 4.27`, but this logarithm shifts by a factor of about
2 with the choice of the inner and outer lengths, so **the agreement claims
nothing beyond the order of magnitude**; `ln(L/l)` itself is not verified.

### 6.2 The end point of spreading against an independently built equilibrium droplet

The spreading is run to equilibrium and the shape it reaches is compared with
an equilibrium droplet **constructed without the solver**, of the same crest
height. Both are measured with the same parabola fit (`--hcut 1.0` µm).

| Nx | dx [µm] | end of spreading | analytic equilibrium droplet | difference |
|---:|---:|---:|---:|---:|
| 500 | 1.2 | 9.4073 deg | 9.4052 deg | **0.02 %** |
| 2000 | 0.3 | 9.4101 deg | 9.4083 deg | **0.02 %** |

The shape itself agrees too: against a crest height of 7.0567 µm the maximum
difference between the profiles is **0.0022 µm = 0.03 %** at `Nx = 2000`
(0.0032 µm at `Nx = 500`). The end point is the steady state of the scheme,
the film of uniform discrete potential ([`numerics.md`](numerics.md) §1.4);
outside the contact line its shape is the continuous equilibrium's even at
`Nx = 500`, where `dx = 1.2` µm is comparable with `ell = 1.15` µm, and the
path there agrees too (apparent angle within 0.06 % of the fine mesh's over the
first second). With the arithmetic mean the disjoining pressure once took for
its face mobility the end point was 1.7 % off at `Nx = 500` and 0.11 % at
`Nx = 2000`. `dx <= ell/5` stays a rule of thumb: this droplet is a lenient
case.

---

## 7. The cost of the coupling

For the block Gauss-Seidel iteration that solves the concentration equation
inside the Newton-Kantorovich iteration of the film-thickness equation,
**the fraction of the iterations that the concentration side governs** was
measured.

### 7.1 Conditions and results

`cases/eq/spread_n500.toml` (Nx = 500, `h* = 0.2` µm, `theta_e = 10` degrees)
with the concentration coupling fully on — convection, diffusion, evaporation
and Marangoni, with a viscosity exponent of 10 — all compared at
`EndTime = 0.1` s. `YUZEN_DEBUG_RESIDUAL=1` prints the (thickness correction,
concentration correction) of every iteration, and the iterations on which the
concentration was the larger are counted.

| Condition | Time steps | Mean NK | Max NK | Conc.-dominated iterations | Dominant at final iteration |
|---|---:|---:|---:|---:|---:|
| Evap. 1e-7 / no disjoining pressure | 1129 | 2.56 | 3 | 57.7 % | 91.5 % |
| Evap. 1e-7 / with disjoining pressure | 1129 | 2.56 | 3 | 60.7 % | 88.0 % |
| Evap. 1e-6 / with disjoining pressure | 1129 | 2.69 | 3 | 96.3 % | 92.3 % |
| Evap. 1e-5 / with disjoining pressure | 1129 | 3.04 | 4 | 98.7 % | 96.4 % |

At evaporation 1e-5 the concentration moves from a maximum of 0.500 to 0.619
(mean 0.500 to 0.571) within the 0.1 s. The stagnation cut-off warned **0
times** in all 4. The corrections at the 1000th step of the two runs at
evaporation 1e-7:

```
  with disjoining pressure     thickness     2.62e-06 -> 7.23e-11 -> 1.84e-15
                               concentration 3.21e-06 -> 9.18e-10 -> 3.51e-14

  without disjoining pressure  thickness     6.20e-06 -> 1.38e-10 -> 3.47e-15
                               concentration 3.70e-06 -> 1.43e-09 -> 7.94e-14
```

The thickness falls by four to five decades an iteration and is at rounding
after two, with the disjoining pressure on or off: `Pihat ~ H^-n` is stiff, but
its Jacobian is exact (§11), so it costs no extra iterations (found taking four
to five before its Jacobian was corrected). **"The fraction of iterations
dominated by concentration" measures not that the concentration is hard but
that the thickness is easy**: the row at 98.7 % takes 3.04 iterations on
average, against 2.56 for the rows with the lowest figures. (Re-measured: the
better first guess of a step that the parallel-efficiency algorithms brought
cut the first thickness correction by two decades, and the mean count with
it; one face mobility for the disjoining pressure moved the percentages by
less than one point.)

The same measurement on the cases that strengthen the coupling
(`cases/stress_*.toml`, all derived from `allphys_1d_reg`); the stagnation
cut-off warned zero times in all four.

| Case | What is changed | Time steps | Mean NK | Max NK | Conc.-dominated iterations |
|---|---|---:|---:|---:|---:|
| `stress_evap` | Evaporation coefficient ×100 | 675 | 3.01 | 4 | 68.0 % |
| `stress_coupling` | Viscosity exponent 10→30; concentration coefficient of surface tension ×10 | 675 | 2.19 | 4 | 2.2 % |
| `stress_dt` | dt ×100, fixed (step control off) | 21 | 4.29 | 6 | 14.4 % |
| `stress_all` | All of the above | 21 | 5.05 | 6 | 0.0 % |

### 7.2 What this measurement does and does not support

**What it supports.** Over the eight conditions the outer iteration converged
in **at most six iterations** and the stagnation cut-off never fired. Because
the convergence test looks at both corrections, the state it converges to is
**the fully implicit coupled solution of thickness and concentration**; writing
the off-diagonal blocks of the Jacobian would change only the number of outer
iterations, which are 2 - 6.

**What it does not support.** It is a measurement of iteration counts, not of
wall-clock time; it is 1-D and `WM` only, and says nothing about 2-D or `CZ`;
and the step counts are those the adaptive step chose (on the spreading family
the step grew without a cut under all four conditions, hence the same 1129
steps), the end time of 0.1 s being the only thing held fixed.

### 7.3 The coupling passes of `CZ`

With the Marangoni term on, `CZ` solves the film and the concentration of a step again at the
midpoint concentration $`c_m`$ the last pass found, Aitken-relaxed, until the midpoint settles (at
most eight passes).

**The gain, derived.** Linearized about a flat film of thickness $`h`$ at rest with a profile
$`c_0(\zeta)`$, a wave $`e^{ikx}`$ sees $`c_m`$ only through the Marangoni stress of its surface
value, and the concentration sees the film only through the $`\zeta`$ velocity that carries
$`\partial c_0/\partial\zeta`$; a plain pass is rank one per wave, and its gain is, for a linear
profile $`\partial c_0/\partial\zeta = \Delta c`$,

```math
\lambda_k = \frac{\Delta t}{2} D_M k^2 \psi(1),\qquad D_M = \frac{h\,\sigma_c\,\Delta c}{2\mu},
```

with $`(1 + \Delta t D k^2/2)\psi - (\Delta t D/2h^2)\psi'' = \zeta(1-\zeta)\varphi`$,
$`\psi'(0) = \psi'(1) = 0`$, $`\varphi = m_k - \frac{a}{1+a}(1 - \zeta/2)`$,
$`a = \Delta t\,\sigma h^3 k^4/(6\mu)`$, $`k^2 = (2 - 2\cos k\Delta x)/\Delta x^2`$ and $`m_k`$ the
velocity's fourth-order surface gradient over the film's two-point one (1 to 7/6). The sign is that
of $`\sigma_c\,\partial c/\partial\zeta`$; with diffusion across the film within the step it tends to
$`\tfrac12\mathrm{Ma}_h\bar s\,\theta/(1+\theta)`$, $`\theta = \Delta t D k^2/2`$,
$`\bar s \in [1/24, 1/6]`$, $`\mathrm{Ma}_h = \sigma_c\Delta c\,h/(\mu D)`$, so it never passes
$`|\mathrm{Ma}_h|/12`$ (×7/6 on the shortest waves). And
$`1 - \lambda_k = \det(I - \tfrac12\Delta t A_k)/[(1+a)\det B_k]`$: the gain passes 1 only where a wave
of the model grows faster than 2/dt. Aitken's factor tends to $`1/(1-\lambda)`$ of the wave that holds
the residual (negative for $`\lambda > 1`$), so it removes one wave, or a band of about one gain,
whatever its gain, but not waves with gains on both sides of 1 together. `tools/cz_coupling.py`
evaluates $`\lambda_k`$ for a grid; `YUZEN_DEBUG_COUPLING` prints every pass of a run.

**Measured.** An evaporating flat film (8.8 µm, Δc = 0.04, D = 1e-9, 64 cells over 0.5 mm, Nz = 11)
whose last step took ten plain passes: the ratio of successive passes against the largest
$`|\lambda_k|`$ predicted for the run's own column —

| σ_c [N/m] | dt = 2 ms | 5 ms | 10 ms | 20 ms |
|---|---|---|---|---|
| +0.03 | 0.336 / 0.358 | 0.836 / 0.843 | 1.419 / 1.529 | 1.955 / 2.636 |
| −0.03 | 0.336 / — | 0.836 / — | 1.419 / 1.539 | 1.955 / 2.602 |
| +0.05 | 0.560 / 0.542 | 1.393 / 1.392 | 2.365 / 2.589 | 3.250 / 4.322 |

(— : the residual reached rounding first; at 20 ms the ratio was still rising). Within 0.97–1.09 for
gains of 0.34 to 2.4. At 10 ms Aitken's factor read −0.64 and 0.28 for σ_c = ±0.05, i.e.
λ = 2.56 and −2.57 (±2.37 predicted), and the relaxed passes took the residual down by 1e4 and 1e6
in nine passes where the plain ones grew it by 2.6 a pass.

On `cz_strat_1d` (fixed step, 0.02 s) the most passes a step took were the same on 32, 64 and 128
cells: 1 at 0.5 and 1 ms, 2 at 2 and 5 ms, 4–5 at 10 ms, not settled at 20 ms (the midpoint still
moved by 0.56–0.68 of the first pass's change); likewise with σ_c of the opposite sign and a fifth of
the diffusivity. With `Fluid.Viscosity.Exponent = 0` every step (32 cells) settled in 2 at 20 ms: what fails
there is the coupling of every level through μ(c) under that state's strong flows, measured, not
derived (a 20 ms step moves its wave by 1.2 rad). With 4× σ_c (a film unstable in the model) up to 8,
not settled from 10 ms. An evaporating flat film under step control (σ_c = ±0.01 to ±0.05,
`Criterion` 1e-2 and 1e-1) took 1–3 passes, up to 8 for σ_c > 0, and one step in 55 did not settle
in two runs; `concz_1d_reg` with evaporation up to 3e-5, σ_c up to ±0.03 and 120–480 cells never
more than 2.

**Found:** a step whose passes did not settle was kept with a warning, unseen by the time-step
control (its truncation-error test looks at the film before the passes): on `cz_strat_1d` at 32 cells
with `InitialDt = 0.02` and a `Criterion` of 1 the error was 2.1 of the amplitude. **Corrected:** such
a step fails like a linear solve that does not converge — cut and solved again under step control
(the first reported in `Warnings.txt`, all counted in `JobInfo.txt` as `Steps failed: the passes did
not settle`; that run ends at 0.17, the error of the step the truncation-error test picks at 0.3),
and a fixed-step run stops naming the passes. A midpoint that moves by less than ten times
`Numerics.ConcEquation.LinearSolver.Epsilon` (and no less than 1e-10) counts as settled. `features`
`coupling_cz` checks both.

---

## 8. The same problem posed twice

Most terms have no known answer, but every term has changes to the input that
must not change the physics. `tools/equivalence.py` (`ctest`: `equivalence`)
runs a shipped case and a copy rewritten in one of those ways, at a fixed time
step, and compares the two point by point; **33 checks**, limit 1e-7, worst
difference relative to how far the field moved:

* **Symmetries, well-mixed**: x ↔ y on a grid with Δx ≠ Δy (2e-9), x → −x
  (1e-15), the same under a strong rotation with UTOPIA for the film height and
  with a fixed off-center nozzle (also on cells 5.6 times longer one way;
  ≤ 2e-9); 2-D uniform and periodic in y against 1-D, also with the Koren
  concentration scheme (0: bit for bit) and with a `FreshInlet` and an `Outlet`
  on either axis and either side (≤ 9e-10); a mirrored 1-D inlet run; a periodic
  film moved 13 cells (1-D) or 5 and 7 cells (2-D) (≤ 2e-9); other scales for the
  dimensionless form (`TotalHeight` × 2.5 and L0 = 0.7 mm, 6e-9); a mirrored 1-D
  nozzle (4e-12).
* **Symmetries, `CZ`**: mirror, shift by 11 cells, transpose, 1-D against 2-D,
  other `TotalHeight` and L0 (≤ 2e-13); a drying 2-D film with the Marangoni
  term transposed and mirrored (3e-11); a column next to a wall restarted at
  c = −0.05 and filled from the columns beside it, transposed (2e-15; it found
  nothing, but a fill that read its neighbors in y from the wrong side passed
  every other check).
* **Against a closed form, the other model or another run**: a uniform `CZ`
  concentration stays uniform while the film moves (0); a drying flat `CZ` film
  keeps its solute (0; over relief 5e-15); the in-plane gravity and centrifugal
  velocities on a flat film are $`\rho g h^2/3\mu`$ and $`\rho\omega^2 r h^2/3\mu`$
  (2e-14); `CZ` carries the solute as `WM` does (0.3 % away from the walls,
  allowed 3 %), which fixes the direction of its advection; and a step rejected
  18 times and retried agrees with the same step taken at its Δt from the start
  (0).

The film thickness itself agrees to 2e-13 or better in every well-mixed check
and in the `CZ` ones without Marangoni flow. The checks found that a retried
step built the explicit half of its Crank-Nicolson update from the rejected
attempt (3.5e-5 of the change), and that the nozzle footprint was measured at
the end of each sub-division of a cell instead of its center; both corrected.

`tools/concz_vs_wm.py` (`concz_vs_wm`, `concz_vs_wm_2d`) poses a film to the
two concentration models with `Fluid.Viscosity.Exponent = 0` and evaporation
and solutocapillarity off, so the concentration cannot affect $`h`$: the peak
thickness agrees to 1e-6 relative and the `Volume` of `History.dat` exactly, in
1-D and 2-D. It found that the two models wrote the volume of a 1-D film for
different strip widths; corrected, both write a strip 1 m wide.

---

## 9. Problems with a closed-form answer

`tools/analytic.py` (`ctest`: `analytic`) runs problems the lubrication
equations answer in closed form. The linear ones follow the Fourier coefficient
of a sinusoid of amplitude 1e-3 of a flat 5 µm film on a periodic domain one
wavelength (0.25 mm) long. Selected results (the manual lists every check):

| Check | What is compared | Result |
|---|---|---|
| `leveling_1d`, `_2d`, `_dt` | capillary decay rate $`\sigma h_0^3 k^4/3\mu`$ | order 1.99 (1-D) and 1.98 (2-D) in Δx, 2.00 in Δt |
| `leveling_cz` | `CZ` against `WM` | 1e-11 |
| `leveling_slip_1d`, `_dt` | Navier slip $`\lambda = h_0/2`$: decay rate $`\sigma(h_0^3/3 + \lambda h_0^2)k^4/\mu`$ | 6.4e-3 at 32 cells; order 1.99 in Δx, 2.00 in Δt |
| `energy_leveling`, `_slip`, `_box`, `_2d` | `Output.Energy` on a 30 % sinusoid leveling (periodic; slip; Symmetry walls with vertical gravity; 2-D): balance of the run $`\vert\mathcal F(T)-\mathcal F(0)+\int\mathcal D\,dt\vert / \int\mathcal D\,dt`$; F never rises from step to step | order 2.00 in Δt (all four; 2e-4 → 5e-5 leveling); rise 0 |
| `energy_dewetting` | a film 5.5 h* thick ruptured onto its precursor (10 s) and gathered into a drop (40 s), 128 cells, Δt = 0.125, 0.0625 s | order 1.97 in Δt; F never rises |
| `energy_saturate_mono` | the film of `spinodal_1d_unstable`, a 30 % sinusoid grown into a steady film (40 s): F never rises, D ≥ 0 | 0 (5.9e-5 with the mean face mobility of the disjoining pressure) |
| `eq_slip` | the equilibrium droplet of §5.3 with slip 1 µm, 0.01 s: max \|Δh\| / max h | 9.5e-6 (6.5e-6 without slip) |
| `gravity_1d` | decay rate with vertical gravity | 6e-3 |
| `spinodal_1d_unstable`, `_dt` | disjoining pressure, $`h_0 = 2h_*`$, a 1 mm sinusoid: growth rate $`\omega = -(h_0^3/3\mu)k^2(\sigma k^2 - \Pi'(h_0))`$ | 2.3e-3 at 32 cells; order 2.00 in Δx and in Δt |
| `potential_1d` | `WritePotential`: the `PHI` column against $`\sigma\nabla^2 h + \Pi(h) - \rho g h`$ on a sinusoid | 1.5e-15 |
| `leveling_logout`, `_times` | the shipped `leveling_1d_logout`, 141 outputs from 1e-7 to 1 s landed on (`ExactTimes`): output time against listed; decay rate at each output | 0 (bit-exact); 1.6e-3, as at a fixed step on 64 cells |
| `spinodal_1d_stable`, `_dt` | the same at $`h_0 = 1.25h_*`$ (below $`1.5h_*`$): decay rate | 3.9e-3 at 32 cells; order 1.99 in Δx, 2.01 in Δt |
| `wave_1d`, `wave_1d_koren` | damping added by the upwinded in-plane gravity | order 3.25 (UTOPIA), 3.41 (Koren) |
| `wave_1d_speed`, `wave_1d_koren_speed` | wave speed $`\rho g h_0^2/\mu`$, 32 cells | 5e-5 (UTOPIA); 6.7e-4, order 2.85 (Koren) |
| `wave_cz` | `CZ` decay rate with in-plane gravity | order 2.10 in Δx |
| `spin_1d` | $`h_0/\sqrt{1 + 4\rho\omega^2 h_0^2 t/3\mu}`$ on a strip | 4e-7 |
| `spin_balance` | a spinning closed 1-D box, the sideways outflow counted | ≤ 1e-12 |
| `evap_wm`, `evap_cz` | drying flat film against the ODE | order 2.00 and 2.16 in Δt |
| `slab_cz` | `CZ` profile against the slab series, Nz = 20 | 3e-4 (second order in Nz) |
| `dry_cz` | `CZ` drying at Eh/D ≈ 1 with D(c), against a reference | order 1.81 in Nz |
| `evap_cap`, `evap_solvent`, `evap_threshold` | the cap at h/2 and the solvent limit $`E\Delta t \le h(1-c)`$ | 3e-8, 3e-8, 1e-7 |
| `dried_wm`, `dried_cz` | a film dried to c = 1: solute kept, c ≤ 1 | ≤ 1e-12 |
| `conc_overshoot`, `_koren` | a front c = 0.99 \| 0.3 carried by the flow, UTOPIA and Central4 (1-D), UTOPIA (2-D), Koren at Courant 1.5: the overshoot past c = 1 moved, not cut | 2e-14; 1e-15 |
| `conc_drying` | the front c = 0.9 \| 0.3 on a closed film drying by Raoult's law as it is carried, UTOPIA and Koren (1-D), UTOPIA (2-D): `ResinVolErr` | 2e-15 |
| `conc_uniform` | a uniform c in a closed box under a nozzle, a strong in-plane gravity and (1-D) the spin: max \|c − c0\| | 6e-16 |
| `conc_undershoot`, `_cz` | a front c = 0.7 \| 0.005 carried by the flow, UTOPIA, KK, Central4 (1-D), UTOPIA (2-D), `CZ`: the undershoot below 0 filled, c ≥ 0, `ResinVolErr` | 1e-14 |
| `wall_balance`, `_2d`, `_pinned`, `_cz`, `_slip` | closed boxes, every term toward the walls, four schemes, every edge and corner (`_slip`: with Navier slip of 3 µm) | ≤ 1e-12 |
| `leveling_adaptive`, `_cz` | adaptive step: decay at the reported times; clock against the sum of the steps solved | 5e-4; 1e-12 |
| `inlet_outlet`, `solute_balance`, `_2d` | a bump carried out through an `Outlet`; the film drying (UTOPIA, Koren, adaptive step, no concentration); 2-D, uniform in y and open on all four sides: max \|`VolErr`\|, \|`ResinVolErr`\| | ≤ 1e-12 (1e-15 to 1.1e-14) |
| `nozzle_balance`, `nozzle_solute` | nozzles 1-D and 2-D, and `nozzle_1d_rotating` thickened a hundredfold | ≤ 1e-12 (1.3e-13) |
| `nozzle_supplied` | the same with the liquid supplied at `Conc` (`ConcType = "supplied"`): `ResinVolErr` counts $`\int c_S S`$, and `VolErr`; on the thick film steps whose Newton correction stagnated are cut | ≤ 1e-12 (2e-13) |
| `nozzle_mixing` | liquid at `Conc` on a film that does not flow, whole domain and a small nozzle, 1-D and 2-D: $`c = (h_0 c_0 + S t\, c_S)/(h_0 + S t)`$; `"local"` keeps $`c_0`$ | ≤ 1e-12 (4.2e-15, rounding) |
| `nozzle_arm`, `nozzle_area` | where a rotating nozzle puts the liquid; what a 2-D nozzle supplies on long thin cells | 2e-4 of the arm; 3.6e-3 |

Found and corrected through these checks: the well-mixed concentration was
enriched at the uncapped evaporation rate (30 % too much solute in one step at
the cap); a drying film could pass c = 1, and the cap then deleted solute (corrected
by limiting the evaporation to the solvent there is); the blend of the `negPow`
profile was not held constant beyond `NegPowRmax`; `History.dat` left out the
solute diffused through the edges and, in 2-D, took the film in the boundary
cell rather than at the face; it then counted what crossed an open edge from the
velocities and face values of the field after each step, not from the face
fluxes the step carried (`VolErr`, `ResinVolErr` 6e-7 on 100 cells, order 1.8
in Δx; corrected: 1e-15); the nozzle arm's ramp times were multiplied by
T0 instead of divided; a 2-D nozzle on cells more than 3.9 times longer than
wide left out cells it partly covered (81 % of its rate on cells 11 times
longer); and in `nozzle_1d_rotating` `VolErr` was −5.4e-5, because the
capillary $`A\cdot\varphi`$ was formed as a product of the entries and the 1-D
`Symmetry` ghosts came off their mirror by rounding. Corrected by forming the
capillary and vertical-gravity parts face by face and satisfying the boundary rows
exactly ([`numerics.md`](numerics.md) §2.1); `VolErr` is 2e-13.

---

## 10. Manufactured solutions

With `Numerics.ManufacturedSolution` the solver adds the source that makes
$`h = h_0 + a\cos(k_x x - \Omega t)\cos(k_y y)`$ an exact solution, with
$`a = 0.3\,h_0`$ (the mobility varies sixfold over a wavelength);
`tools/mms.py` (`ctest`: `mms`) requires order ≥ 1.8 in Δx and in Δt. The
sources are derived with SymPy by `tools/mms_generate.py` from the equations as
the manual states them (`mms_source_current` checks that the generated file is
current); the evaporation is written out by hand.

| Case | In Δx | In Δt |
|---|---|---|
| `WM`, 1-D / 2-D | 2.61 / 2.78 | 2.08 / 2.15 |
| `CZ`, 1-D / 2-D | 2.02 / 2.01 | 2.06 / 1.94 |
| `WM`, 1-D / 2-D, every term | 1.89 / 2.11 | 2.07 / 2.07 |
| `CZ`, 1-D / 2-D, every term | 2.00 / 2.05 | 2.05 / 2.01 |
| `CZ`, c(x) / c(x,y) | 2.00 / 1.97 | 2.22 / 2.32 |
| `CZ`, c(x,ζ) / c(x,y,ζ), stratified | 2.01 / 1.96 | 2.30 / 2.02 |
| `WM`, 1-D, Navier slip / every term with slip | 2.90 / 2.03 | 2.12 / 2.07 |

"Every term" means capillary pressure, both gravities and the centrifugal term
about an off-center axis over a 2 µm relief, and in the well-mixed model also
the disjoining pressure, the Marangoni flux, $`\mu(c)`$, $`\sigma(c)`$, $`D(c)`$
and the evaporation (`PowerN` law in 1-D, Raoult with the `negPow` profile in
2-D), with a manufactured concentration. The well-mixed orders above two are
the pre-asymptotic mixture of second-order capillary and third-order upwinded
terms. The stratified `CZ` cases add $`0.1\cos(\cdot)\cos\pi\zeta`$ through the
film and check $`h`$, $`c`$ and $`\mathbf{u}`$ at every level; evaporation and the
disjoining pressure stay out of `CZ`, and the nozzle out of both models.

The slip rows (`wm_slip_1d`, `wm_all_slip_1d`) take a slip length of half the
mean film, so that $`\lambda h^2`$ is as large as $`h^3/3`$; their source is
generated for the slip mobilities. `wall_balance_slip` found that the face
mobility of the vertical gravity has to vanish against the ghost $`-h`$ of a
`Pinning` edge, as the mean of $`h^3`$ did by symmetry and the mean of the slip
mobility did not (VolErr 2.5e-6); the vertical gravity now takes the positive
face mobility, which is 0 there with and without slip.

These runs found, all corrected: the body forces and the vertical gravity over
relief took the viscosity of the cell rather than of the face (not
conservative); the well-mixed concentration formed its first explicit half
before the first velocities (first order in Δt); `CZ` did not count the 1-D
centrifugal outflow; in `CZ` the capillary term used one face's surface tension
for both, the centrifugal x velocity ignored the axis offset, the first Simpson
interval was misplaced, the 2-D capillary velocity lacked its cross
derivatives, and seam faces took the far side's velocity (as did the
well-mixed model); and with stratification, the film-height flux took $`\sigma`$
and $`\mu`$ at the depth mean, `Uz` was stored one cell high and was not the
vertical velocity, and the surface ghost was not the image of the cell below.
A step limited as Δx² by two disagreeing velocities was removed by forming the
velocity from the film of the step itself and iterating the coupling to the
midpoint: an evaporating 1-D film that needed 1548 adaptive steps needs 134.

---

## 11. Jacobians and the convection schemes

`test_jacobian.cpp` compares, term by term in 1-D and 2-D, $`Jv`$ with a central
difference of $`R(h) = A(h)h`$ on a non-flat film over non-flat relief with a
varying concentration (rounding 1e-10; allowed 1e-7), with seven convection
schemes for the three upwinded terms. It found three faults, corrected: the
vertical gravity and the disjoining pressure took the derivative of the face
mobility at the neighbor for the cell's own column (2e-2); the disjoining
pressure added all of $`\hat\Pi'`$ (50 % off); and the 1-D centrifugal outflow
added $`2\,\mathrm{Sp}\,h/\mu`$ where $`2\,\mathrm{Sp}\,h^2/\mu`$ was due (3e-2).
The same file checks the face fluxes the well-mixed concentration is carried with
(`FaceFlux`): every term and scheme, 1-D and 2-D, their signed sum per row is the
row's $`A h`$ to 1e-13, and the two cells of a face see the same flux.

`test_finite_differences.cpp` fixes, for every scheme, exactness on
polynomials, order, the bounds of the limited schemes, mirror and sign
symmetry, the signs of the frozen weights and the derivative.
`tools/advection.py` (`ctest`: `advection`) carries a step and a sine across a
periodic flat film at a Courant number of 0.4, in both models:

| Scheme | Overshoot of the step, WM / CZ | Sine, maximum error at 128 cells (order) |
|---|---|---|
| KK | 6.3 % / 7.0 % | 1.9e-4 (2.99) |
| UTOPIA | 8.6 % / 10.4 % | 6.3e-5 (2.99) |
| Central4 | 36 % / 39 % | 6.4e-4 (2.00) |
| Koren | 0 | 5.1e-3 (1.56; mean error 2.35) |
| van Leer | 0 | 1.0e-2 (1.41; mean error 1.98) |

The bound holds for a Courant number up to one only (Koren overshoots 16 % at
1.5); `advection.py` also checks that the solver's Courant report fires where it
should.

The sine is carried along x without diffusion. `oblique_wm` (part of
`advection`; alone with `--oblique`) checks the order of the well-mixed
concentration in 2-D with convection and diffusion together: a flat square
periodic film, 1 mm across, carried at 1 mm/s in the direction (2, 1)/√5, with
$`c = 0.5 + 0.2\sin k(x+y)`$ oblique to the flow and to both axes and
$`D = 5\times10^{-7}`$ m²/s, on 16, 32 and 64 cells a side with 20, 40 and 80
steps over one decay time (Δt refined with Δx). From the modified equation the
error is second order (the central diffusion under-damps by
$`D(k_x^4+k_y^4)\Delta x^2/12`$) with a third-order part of the other sign
(UTOPIA over-damps by $`(|u_x|k_x^4+|u_y|k_y^4)\Delta x^3/12`$), whose share
$`r`$ is 0.67 times the cell Péclet number; the observed order is
$`2 + \log_2[(1-r)/(1-r/2)]`$, below 2 until $`r`$ is small. This cancellation
is what keeps the order of $`c`$ low in the 2-D well-mixed manufactured solution
with every term on at the 1-D diffusivity (§10). At
cell Péclet numbers of 0.125 to 0.031 the formula gives 1.94 and 1.97; UTOPIA
converges at 1.91 and 1.97 (error within 3 % of the formula's on every grid)
and Koren at 1.88 and 1.92. Both must reach 1.9–2.1 on 32→64 cells, with the
order on 16→32 within 0.1 of it. First-order upwinding gives 0.72 and 0.87,
KK 1.76 and 1.90. The film stays flat, so the film height is not checked here.
The six runs take about 20 s on two cores.

On a steep block sliding over a 0.1 µm precursor the lowest film height
is 0.085 µm with KK, 0.028 µm with UTOPIA and 0.0998 µm with Koren and van Leer;
over 0.02 µm only the limited schemes finish. In the manufactured solutions a
limited concentration scheme costs up to a factor of 2 (`CZ`) and 6–7
(well-mixed 2-D) in error against UTOPIA, while Koren for the film height costs
1.2–1.3 times KK's error. Every scheme runs the shipped cases with the same
steps, and with the same Newton iterations except where the step is a hundred
times longer (Koren and van Leer for the concentration: 142 and 139 against
KK's 104 in `stress_all`), and gives the same answer with 1 and 2 threads to
the last bit. Hence the defaults: **`Koren` for the film height, `UTOPIA` for the
concentration** ([`numerics.md`](numerics.md) §1.5). The comparison also found
the Marangoni term upwinded by the sign of Δc rather than of the flux (wrong
with a negative `Fluid.SurfaceTension.Coef`); corrected. A `CZ` film stratified
with $`\sigma_c\,\partial c/\partial\zeta > 0`$ and a weak diffusivity is
unstable in the model at short wavelengths (`tools/cz_stability.py`); the
stratified manufactured solutions use a diffusivity at which every resolved
wave decays.

---

## 12. Restart, features, input guards and error paths

* **`restart_roundtrip`** restarts each model, 1-D and 2-D, from each format and
  compares the first field the restarted run writes with the one it started
  from, with its time and file number. **`restart_chain`**
  (`tools/restart_roundtrip.py --chain`) stops thirteen cases of both models —
  adaptive step, evaporation, the `CZ` coupling passes, a limited scheme, a
  moving nozzle (its liquid at the local concentration, and at its own), probes, the steady-state and dry-out stops, a film dried against a
  `Pinning` wall, a film fed through a `FreshInlet` and drained through
  `Outlet`s (1-D and 2-D) — by `MaxSteps` at steps that are not output times and
  restarts them from the last file, alternating `.vts` and `field*.dat`: every
  field file, every row of `History.dat`, `CalcInfo.dat` and the probes, and
  the step counts must be identical, bit for bit. This rests on the
  `field*.state` file written beside every field file (time to the bit, step,
  the last two steps, the adaptive-step counters, the output schedule, the
  running totals of `History.dat`, the steady-state count, and the fields with
  their ghosts). Without it two fixed-step cases must keep the output times and
  agree to 1e-11 (well-mixed) and 1e-7 (`CZ`). Five faults found and
  corrected: the output schedule restarted from t = 0; the totals started from
  zero; `CZ` took its first step from the initial field; the outer ghost column
  of a `CZ` field read from `field*.dat` kept the initial film; and the last
  step of each part was not counted toward the steady state.
* **`features`** (`tools/features.py`): probes land on the nearest point and
  name it correctly (they did not, and in `CZ` landed nowhere near; corrected);
  the steady-state, dry-out and `MaxSteps` stops; the adaptive step against a
  fixed small step (1e-4 at `Criterion = 1e-4`); `adaptive_clock_*`, the
  reported time against the sum of the steps solved (1e-12, both models, 1-D
  and 2-D); the `CZ` coupling passes driven past their limit of eight (a fixed-step run stops, a
  step-controlled one cuts the step; §7.3); `stagnation_wm`, the first step of `stress_wetting`,
  whose Newton correction stagnates far from convergence (a fixed-step run stops naming the
  stagnation, where taken with a warning the film went negative and the run ended as diverged; a
  step-controlled one cuts the step thirteen times, warns once, counts the cuts and keeps `VolErr` at
  2e-15); `evap_precursor_wm` and `evap_precursor_cz`, a film thinning below
  its precursor while it evaporates, whose `EvaporatedVol` must never fall (the
  precursor correction was negative there and condensed solvent onto the film,
  43 % of the initial volume over 1000 s of `allphys_1d_reg`; corrected); a field
  file edited after it was written (the restart takes it and warns); a restart
  past the `Pinning` check (a warning, not a refusal).
* **`input_guards`**: a misspelled key, a key in the wrong table, a one-sided
  `Periodic`, a `Periodic` edge with the centrifugal term, a `CZ` edge the model
  cannot impose, a pinned wall facing a film of full thickness, a fixed step
  whose Newton or concentration solve fails or whose Newton correction stagnates far from convergence (stops the run, naming the cause)
  and the same under step control (cut and retried), a concentration read
  outside [0, 1] or a film read below zero (stops the run as diverged), a nozzle concentration `Physics.Source.Conc` outside [0, 1]
  with `ConcType = "supplied"` (and an unknown `ConcType`), the grid-size limit $`(N_x+4)(N_y+4)(N_z+4) \le 5\times10^{7}`$,
  and an unterminated `.inp` string are refused with the reason.
* **`error_paths`**: a wrong command line, input file or result directory, a
  missing or truncated restart file, a file size limit or a full device,
  `SIGINT`, `SIGTERM` or `SIGKILL` in the middle of a run, and paths with spaces
  and Japanese characters must stop the run with the path and the reason and
  leave no file cut short. Every field, state and `.pvd` file is written under a
  temporary name and renamed when whole; a failed write stops the run with exit
  code 1; the first `SIGINT` or `SIGTERM` finishes the step, writes the state
  reached and exits with 130 or 143; and the run continued from that state is,
  bit for bit, the run that was not interrupted. Exit code 0 means one of the
  normal stops (or `--check-input`, `--help`, `--version`).
* **`doc_keys`** holds the documentation of every key to the parser (§2).

---

## 13. Fuzzing, random inputs, compilers, sanitizers and long runs

**Fuzzing.** The readers of the input file (TOML and `.inp`) and of the restart
files (`field*.dat`, `.vts`, `field*.state`) were fuzzed with libFuzzer under
ASan and UBSan (about 2.7 h of one core, 14.5 million inputs;
`-DYUZEN_BUILD_FUZZERS=ON`, Clang only). Found and corrected: a toml++
unreachable branch on a non-ASCII bare key; a `.vts` block length trusted
beyond the grid (heap overflow); short ASCII rows read uninitialized; a huge
array length in the state reader; conversions that accepted trailing text,
`nan`, fractional counts, lists for scalars and dates; and values that ran to
nonsense (a negative `MaxDtFactor`, a negative viscosity, a film with no liquid
that ended with exit code 0 and NaN in `History.dat`). All are refused, and the
well-mixed model stops on a field that has broken down, as `CZ` does.

**Random valid inputs.** `tools/proptest.py` draws inputs from a seed (model,
dimension, terms, parameters over orders of magnitude, edges, initial shape,
scheme, fixed or adaptive step) and checks clean ends, finite fields, a
positive film, a bounded concentration for a bounded scheme at Courant ≤ 1,
volume (1e-9) and solute (1e-11 for the well-mixed model with its advection
on, 1e-13 for `CZ`) on closed domains and, well-mixed, through `FreshInlet` and
`Outlet` edges as well, a uniform concentration without
evaporation kept uniform (1e-10; 1e-12 for `CZ`), and, for some,
transpose, mirror, thread and restart-chain invariance. About 990
configurations were run; `ctest` runs 20 small seeds and fifteen regression
seeds. Found and corrected, among others: the well-mixed body forces and
Marangoni stress put liquid through `Symmetry` and `Pinning` walls (`VolErr`
2.5e-4; 1e-14 after the correction), and so did the `CZ` body forces; `ResinVolErr` did not
count the nozzle's solute; a non-flat drying `CZ` film lost 3–25 % of its solute
through the slope part of the diffusive flux, and its concentration equation
was not conservative (corrected by solving for $`hc`$, solute kept to 1e-14); a drying film
passed c = 1; **the time-step control lengthened the step before the clock was
advanced, so the clock ran ahead by up to 0.5 %** (corrected: the longer step applies
from the next step, the clock advancing by the step solved); the 1-D
centrifugal outflow was not counted in the balances; the `Pinning` check was
repeated on a restart; a `Periodic` edge was accepted with the centrifugal
term; and the cap at c = 1 of the well-mixed concentration cut the solute an
unbounded scheme carried past 1 at a sharp front (0.72 % over one crossing of
100 cells with UTOPIA, 3.0 % with Central4), and froze the limiter of Koren and
van Leer at the capped values past Courant 1 (+1 %) — corrected: the excess is
moved to the neighbors with room and the limiter is frozen at the iterate as
solved (`conc_overshoot`, `conc_overshoot_koren`). The well-mixed concentration
equation, solved for c in non-conservative form, kept the solute only to the
accuracy of the scheme (1e-5 to 3e-4 of it on a closed drying front, 50 to 200
cells; proptest allowed 1e-2); corrected: it is the balance of h c with the face
fluxes of the film-height equation ([`numerics.md`](numerics.md) §1.7), solute
kept to 1e-14 over 200 seeds (`conc_drying`, `conc_uniform`). An unbounded
scheme's undershoot below c = 0 at a front was left as it was and could pass
the divergence test of −0.1 (QUICK, Central4; `CZ` UTOPIA); corrected: the
deficit is filled from the cells around holding solute (`conc_undershoot`).
Through a `FreshInlet` or an `Outlet` the balances held only to the error of the
scheme, `History.dat` counting what crossed from the velocities and face values
of the field after each step (`ResinVolErr` up to 1.1e-2, seed 1033); corrected:
the step's own face fluxes are counted, 4e-15 over seeds 1000–1099, and proptest
checks the balances of open domains as of closed ones.

**Compilers.** The shipped cases and `ctest` were run with GCC 13 and Clang 18 at
`-O0`, `-O2` and `-O3`, with and without `-fno-strict-aliasing`, with
uninitialized locals patterned, with `MALLOC_PERTURB_` and with
`-D_GLIBCXX_DEBUG`: every one gives the same field files to the bit, and GCC
and Clang agree. GCC 11 (RHEL 9) builds it and passes `ctest`; it found a header
using `size_t` without `<cstddef>`, corrected. Fused multiply-add and Eigen's vector width change $`h`$ by
about 1e-12 of its maximum (up to 1e-8 on the contact-line cases over
thousands of steps). `-ffast-math` fails three unit tests, and
`-funsafe-math-optimizations` should not be used.

**Sanitizers and coverage.** ASan and UBSan (also in CI on every push), Valgrind
memcheck, TSan with Archer, floating-point traps (`feenableexcept`),
clang-tidy and coverage (90 % of lines) were run in four rounds, and the
sanitizers, Valgrind, TSan and the fuzzers (6.2 million inputs) once more on the
released input format. They found, all
corrected: a retried step built from the rejected attempt; a van Leer limiter
that underflowed for tiny and overflowed for huge differences; tests copying
uninitialized flags; the two rows of an `Outlet` or `FreshInlet` (and of a
1-D pinned concentration) left off by part of the solve's residual; and a
`History.dat` that is a device without end (`/dev/full`) read whole to resume
it, which took all memory before the write failure could be reported. The traps
otherwise fire only where harmless (the initial-shape sigmoid, the 0/0 velocity
of a pinned wall face that the boundary condition overwrites, deliberate
overflow checks, rejected step attempts). No race was found at two to four
threads, and the field files are bit-identical at one, two and four threads.

**Long runs and large grids.** Shipped cases run far longer than they are,
single-threaded:

| Run | Steps | t [s] | max \|VolErr\| | max \|ResinVolErr\| | Memory |
|---|---|---|---|---|---|
| `allphys_1d_reg`, dt = 1e-4 s fixed | 500,001 | 50 | 1.1e-14 | 2.3e-14 | 8.0 MB |
| `allphys_1d_reg`, adaptive | 100,860 | 1000 | 8.6e-15 | 9.7e-15 | 7.9 MB |
| `concz_1d_reg`, dt = 1e-4 s fixed | 200,001 | 20 | 4.0e-15 | 5.3e-15 | 8.1 MB |
| `allphys_2d_reg`, adaptive | 3,612 | 30 | 1.6e-14 | 1.8e-14 | 14.5 MB |

The balances stay at rounding and do not drift. The resident memory is no
larger at the end than a quarter of the way through. A run writing
30,001 field files keeps a valid `.pvd` with distinct time stamps. The peak
memory of a step on larger grids, run for a few steps each:

- well-mixed 2-D at 1024² cells: 3.8 GB with `eigen`, 3.2 GB with `amgcl`;
- well-mixed 1-D at Nx = 10⁶: 2.0 GB;
- `CZ` 2-D at 200 × 200 × 20: 1.6 GB.

That is about 3.7 kB per cell in 2-D, 2.0 kB in 1-D and 2.1 kB per cell in
`CZ`. `python3 tools/long_runs.py --exe build/bin/yuzen --workdir runs/long`
repeats all of them (about an hour on one core).

---

## 14. Mutation testing

To ask whether the checks would notice a fault, single changes (a sign flipped,
an index moved by one, a factor changed, `+=` made `=`, a loop bound moved) were
made to the sources one at a time and `ctest` run on each.

* **Sample 1** (108 changes to the physics sources): 52 caught at first; the
  tests written to close the gaps (the Jacobian test, every-term manufactured
  solutions, nozzle, evaporation-cap, inlet/outlet and drying checks,
  `features`) brought it to 24 more caught. Of the 32 left, none is behavior
  left untested: 13 out-of-bounds accesses (caught by the sanitizer build), 11
  that cannot change a result, 2 that change only the constant of an error, and
  6 in messages and output. They found the `negPow` blend oscillating beyond
  its outer radius.
* **Sample 2** (drawn uniformly, nothing written until it had run: 100 of 3,344
  possible changes in model and numerics, 25 of 312 in input and output):
  68 of 81 effective changes caught (84 %; 89 % with the sanitizer), and 7 of
  25 in the output code. The checks then added brought them to 73 (77) of 81
  and 22 of 22, and found the probe header of a 1-D run naming the wrong point.
* **Sample 3** (after much of the code was rewritten: 100 of 3,815 and 25 of
  309): 74 of 83 caught (89 %; 92 % with the sanitizer) and 20 of 23 (87 %);
  weighted, 89 % and 91 %. The checks then added brought them to 77 (79) of 83
  and 21 of 23, and found two faults, corrected: the 2-D nozzle on long thin
  cells and the unterminated `.inp` string.
* **Sample 4** (the lines changed after sample 3: 30 of 511 possible changes):
  16 of 20 effective changes caught (80 %; 85 % with the sanitizer). The check
  then added, a `CZ` column below c = 0 next to a wall against its transpose,
  brought it to 19 of 20; the one left changes only a debug print. No fault
  found.

---

## 15. Procedures for reproduction

Everything is run from the root of the repository.

### 15.1 The test suite (§2–§5)

```sh
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build -j
ctest --test-dir build                        # all seventeen entries
./build/bin/yuzen_tests                       # case and assertion counts
./build/bin/yuzen_tests -sf="*test_grid.cpp"  # the per-file breakdown (the table in §2)
YUZEN_DBG_EQ=1 ./build/bin/yuzen_tests -tc="*residual velocity on the equilibrium shape*"
```

The last line prints the in-test residual velocities of §5.2 (42.996 / 10.801 /
2.7116).

### 15.2 Thread count and backend (§4.2, §4.3)

```sh
./tools/run_case.sh cases/leveling_1d_reg.toml runs/reg_th1 1
./tools/run_case.sh cases/leveling_1d_reg.toml runs/reg_th2 2
python3 tools/compare.py runs/reg_th1 runs/reg_th2            # => 0.000e+00

./tools/run_case.sh cases/allphys_2d_reg.toml runs/2d_th1 1
./tools/run_case.sh cases/allphys_2d_reg.toml runs/2d_th2 2
python3 tools/compare.py runs/2d_th1 runs/2d_th2
```

The exact agreement of §4.2 is `eigen`'s. On a build with AMGCL, `"auto"`
resolves to `amgcl` for the 2-D `allphys_2d_reg`, and the second pair then
may differ in the low-order bits (§4.3): write `"eigen"` into both `Backend` keys,
or configure with `-DYUZEN_WITH_AMGCL=OFF`. To change the backend, write
`"eigen"` or `"amgcl"` into both `Numerics.HeightEquation.LinearSolver.Backend`
and `Numerics.ConcEquation.LinearSolver.Backend` (no shipped case sets it), and
check the `Solver Info` section of `JobInfo.txt` for what ran. The total
BiCGStab iteration count is the sum of the `BCG_iter` column of
`CalcInfo.dat`.

### 15.3 The equilibrium droplet (§5.2, §5.3, §5.5)

The ASCII restart file is built with `tools/make_equilibrium_1d.py`; the files
in `cases/eq/*.dat` are byte-identical to its output. `--sweep` prints the
table of §5.5 and writes no file, but `--out` is still required.

```sh
python3 tools/make_equilibrium_1d.py --out /dev/null --sweep

python3 tools/make_equilibrium_1d.py --out cases/eq/eq_n2000.dat \
        --nx 2000 --xmin=-0.3e-3 --xmax=0.3e-3 \
        --hstar 0.02 --theta 10 --hmax 1.0

for n in 1000 2000 4000; do
  ./tools/run_case.sh cases/eq/eq_n$n.toml runs/eq_n$n 1     # §5.2
done
```

Dividing the maximum absolute value of the 3rd column (`Ux` [m/s]) of
`field00000.dat` by the velocity scale `U0 = L0/T0`,
`T0 = mu L0^4 / (sigma h0^3)`, gives the dimensionless residual velocity; here
`L0 = 1.749042e-03` m, `T0 = 3.119463e+04` s, `U0 = 5.606869e-08` m/s (also
printed in each run's log under `Scaling factors`).

The time evolution (§5.3) uses `cases/eq/eq_run.toml` with
`SteadyStateCriterion = -1` added to `[Numerics]`; `EndTime = 0.0102` gives the
5 `field*.dat` files at `t = 0, 2.5e-3, ..., 1.0e-2`, and a sixth at the end
time. The apparent angle is measured with `fit_cap` of §15.4.

### 15.4 Spreading (§6)

```sh
./tools/run_case.sh cases/eq/spread_n500.toml runs/spread 2
python3 tools/measure_spreading.py runs/spread/spread_n500 --dt-out 1e-3 --stride 5
```

`--dt-out` matches `Output.FieldData.TimeInterval` in the input. `--stride 5`
uses every 5th file; the correlation is 0.984 with `--stride 1` and 0.999 with
`--stride 10`, so quote the stride with it.

For §6.2, extend `EndTime` to 8.0 and let the steady-state test stop the run
(`t = 2.23` s at `Nx = 500`, `t = 2.24` s at `Nx = 2000`), read the crest
height from the last `field*.dat`, build a droplet of the same crest height
with `make_equilibrium_1d.py --theta 10 --hmax <crest height/h0>`, and put both
through `fit_cap()` of `measure_spreading.py`.

```python
import sys; sys.path.insert(0, "tools")
from measure_spreading import fit_cap
R, theta, hmax = fit_cap("runs/spread/spread_n500/field02232.dat", 1.0)
```

### 15.5 The cost of the coupling (§7)

Copy `cases/eq/spread_n500.toml` and set `Physics.ConcField.EnableConvection` and
`EnableDiffusion`, `Physics.SolutoCapillarity.Enable` and
`Physics.Evaporation.Enable` to `true`, `Physics.Evaporation.Coef` to `1e-7` /
`1e-6` / `1e-5`, `Physics.DisjoiningPressure.Enable` to `true` / `false`, and
`Numerics.EndTime` to `0.1`.

```sh
YUZEN_DEBUG_RESIDUAL=1 ./build/bin/yuzen case.toml 2> resid.txt
```

Split the `[resid] <iteration> <thickness correction> <concentration correction>`
lines of `resid.txt` into time steps at each iteration number 0, and count the
fraction of lines on which concentration > thickness. The `stress_*` cases are
used as they are.

### 15.6 Sanitizers

Configure with `-DCMAKE_CXX_FLAGS="-fsanitize=address,undefined"` and `-DYUZEN_PREFER_LIBOMP=OFF`
(LLVM's `libomp` keeps its thread structures to the end of the process, which LeakSanitizer reports
as leaks) and run `ctest` in that build; run `valgrind --error-exitcode=77 build/bin/yuzen
case.toml` with `Numerics.MaxSteps` set small; build with `clang++
-fsanitize=thread` against `libomp` and run with `OMP_TOOL_LIBRARIES` pointing
at `libarcher.so` (without Archer, TSan reports the OpenMP runtime's own
synchronization as races).

---

## 16. What has not been verified

None of these is a known defect, but each is a place where Yuzen offers no
evidence and you should supply your own.

**Physics and model.**

* **2-D (and axisymmetric) contact lines.** Everything in §5 and §6 is 1-D. The
  first step would be §5.2 repeated on an analytic axisymmetric equilibrium
  shape.
* **Marangoni flow and substrate relief** are verified by the manufactured
  solutions only (§10), not against a closed form.
* **`CZ` is verified against the well-mixed model, not against an analytic
  solution**, apart from its through-thickness diffusion (slab series), its
  advection (against `WM`) and the drying profile of a flat film (`dry_cz`).
  Such a profile on a film that is not flat is not checked.
* **A restart without a state file is an approximation**: the adaptive step,
  the totals of `History.dat` and the steady-state count start afresh.
* **Contact-angle hysteresis** and the regime in which the dynamic angle departs
  from Cox-Voinov are not treated.

**Numerics.**

* **The order in time is measured with every term on, not for every
  combination of switches.** In 2-D the well-mixed manufactured solution is
  checked at a tenth of the 1-D diffusivity; the asymptotic second order of
  the 2-D concentration is checked on an oblique sine carried and diffused on
  a flat film (`oblique_wm`, §11), and the approach to it in the manufactured
  solution at the 1-D parameters was measured once, outside `ctest`.
* **The limited schemes are bounded only for a Courant number up to one.** The
  solver reports a step past it but does not shorten it; `nozzle_1d_rotating`
  (4.6) passes it in the film-height equation, and `stress_all` reaches 0.94. Newton's
  convergence at Koren's kinks was checked on the cases of §11 only.
* **The convergence of the `CZ` Marangoni coupling is derived for a flat film
  only** (§7.3: the gain of a pass, measured to within 10 %, passes 1 only where
  a wave of the model grows faster than 2/dt); the coupling through a viscosity
  that depends on c is measured on one manufactured state, not derived. A step
  that does not settle in eight fails and is cut (or stops a fixed-step run);
  how often in long runs other than those measured is not known. A stratified `CZ` film with a weak diffusivity is
  unstable in the model at short wavelengths, and no check says what a
  converged answer there would be. Whether the truncation-error estimate of
  the adaptive step corresponds to the actual error is checked only by
  `adaptive_dt`.
* **The `amgcl` backend is outside the test suite** apart from passing `ctest`
  on an AMGCL build; its agreement and thread differences (§4.3) come from one
  case run by hand, and it is not bit-reproducible by construction.
* **`FreshInlet` and `Outlet` are checked against each other** (mirrors,
  transposes, 1-D against 2-D, balances), not against a known inflow.

**Performance and scale.**

* Wall-clock time and parallel scaling are not verification measurements; the
  timings in [`numerics.md`](numerics.md) §7 are indicative.
* The meshes used reach `Nx = 4000` in 1-D and `Nx = Ny = 120` in 2-D. Larger
  grids were run for a few steps to measure their memory (§13); how the linear
  solver converges on them and how rounding accumulates over a long run on
  them have not been measured.
