# The example and regression cases

**The cases are the `.toml` files.** That is what to read, copy and run.

```sh
mkdir -p run && cd run && cp ../cases/tutorial_1d.toml .
../build/bin/yuzen tutorial_1d.toml
```

All cases are TOML, which is the input format the manual specifies. Every case
sets `JobName` to its own file name, so a run writes its results into a
directory named after the case.

> The solver also reads a block-structured `.inp` format — it picks the parser
> from the extension — which the manual does not specify.

The manual's *The example cases* appendix has the full table, including
single-threaded run times and the warnings each case is expected to raise.

## The 1-D cases — regression, tutorial and nozzles

1 mm domain, `Nx = 300`, a 10 µm height scale, an elliptical bump on an
elliptical substrate feature, well-mixed concentration.

| | |
|---|---|
| `leveling_1d_reg` | Capillary leveling, Marangoni, evaporation. The baseline. |
| `allphys_1d_reg` | The same with every term on, including gravity and the centrifugal term. **Run this first after changing the solver.** |
| `allphys_1d_neg_reg` | As above with the body-force signs reversed, so a sign error cannot cancel. |
| `periodic_1d_reg` | Periodic in *x* instead of symmetric, the step made a centered plateau so that the ends meet. |
| `leveling_1d_logout` | A small sinusoid levelling on a periodic film, no concentration, the field written at 141 times spaced logarithmically over seven decades (`Output.FieldData.LogSchedule`). |
| `concz_1d_reg` | The `CZ` model, `Nx = 60`, concentration resolved through the thickness. |
| `tutorial_1d` | The tutorial case in the manual's *Your first calculation*. |
| `smoke_1d` | The quickest case here. A smoke test. |
| `nozzle_1d_fixed` / `nozzle_1d_rotating` / `nozzle_1d_offset` | Liquid supplied through a nozzle: stationary arm, rotating arm, nozzle off center. |

## The 2-D cases — the same physics in plan view

These are what would catch a mistake in the 13-point stencil, in the two
direction orderings, or in the row decomposition. Their initial condition is a
straight front crossing the domain at about 20° to the *x* axis — the domain
sits on the flank of an ellipse far larger than itself — so nothing about the
answer is aligned with either grid direction and *y* carries as much variation
as *x*. With a feature symmetric about *y*, a 2-D case would reproduce a
1-D answer, and an error in *y* could go unseen.

The front is softened over about 100 µm: six cells at 60², four at 40² and
twelve at 120². `SlopeWidth` is what sets that, and for an `EllipseBump` it is
measured in the **normalized radius**, not in meters — write it as a length and
the front lands inside a cell or two. The manual's `EllipseBump` section gives
the formula.

| | |
|---|---|
| `leveling_2d_reg` | 60², the 2-D regression case. |
| `leveling_2d_tight_reg` | 60² with tighter tolerances — the answer should not move. |
| `leveling_2d_fine_reg` | 120², the largest shipped and the slowest to run. |
| `allphys_2d_reg` | 40², every term on, coarse so that it stays quick. |
| `concz_2d` | 60²×30, the `CZ` model on a two-dimensional film. Not a regression case — it is the worked example of the pseudo-3D combination, and takes about 70 s single-threaded. |

## `stress_*` — hard coupling

All four derive from `allphys_1d_reg` and exist to push the
height–concentration coupling. None of them raises a stagnation warning.

| | |
|---|---|
| `stress_evap` | Evaporation coefficient ×100. |
| `stress_coupling` | Viscosity exponent 10 → 30, surface-tension coefficient ×10. |
| `stress_dt` | Time step ×100, adaptive stepping off. |
| `stress_all` | All three at once. |

`stress_wetting` is the other kind: a slab meeting a thick precursor
(`h*/h0 = 0.2`) on a coarse mesh, so that the first steps stagnate and are
cut. It tests the Newton iteration and the step control, and is not a model
of a contact line.

## `wet_*` — contact lines

Disjoining pressure on, so there is a real contact line.

| | |
|---|---|
| `wet_static` / `wet_static_5` / `wet_static_15` | A 20° droplet spreading to equilibrium at θe = 10°, 5°, 15°. |
| `wet_coupled` | A contact line with the full concentration coupling. |

All four use a precursor film of `h*/h0 = 0.01`, resolved by a mesh of
0.2 µm. A droplet on a precursor that thin cannot start from a slab: where
the slab meets the precursor the velocity is many orders of magnitude above
the natural scale and the step collapses. So each starts from an
*analytically constructed equilibrium droplet* supplied through the restart
mechanism, built by `tools/make_equilibrium_1d.py` with the command recorded
in the case: `wet_static_drop.dat` (20°) and `wet_coupled_eq.dat` (10°).
**The `.dat` files are part of the cases, not outputs** — copy them
alongside the input or the restart fails.

## `eq/` — the verification cases

The inputs behind the disjoining-pressure and Cox–Voinov sections of
`docs/validation.md`. All start from an equilibrium droplet built by
`tools/make_equilibrium_1d.py`, and the `.dat` files here are byte-for-byte
that script's output. They use `h*/h0 = 0.02`, so that `Hmax/H* = 50` and the
apparent angle is within 1 % of the requested one.

| | |
|---|---|
| `eq/eq_n1000` / `eq_n2000` / `eq_n4000` | Residual velocity at Δx = 0.6, 0.3, 0.15 µm. Together they give the measured order 1.98. |
| `eq/eq_run` | Time evolution from the equilibrium shape. Switch the steady-state test off or it stops at once. |
| `eq/spread_n500` | Spreading from a 20° droplet at `ContactAngle = 10.0`; the Cox–Voinov measurement. |

## `scaling/` — the SMP scaling benchmark

Generated, not hand-written: `tools/scaling/make_cases.py` writes them together
with `manifest.tsv`. **Do not edit them by hand** — edit the table at the top
of the generator and regenerate. They differ from every case above in fixing
the work per run (constant time step, fixed step count, no field output) so
that wall-clock times are comparable between thread counts. See
`tools/scaling/README.md`.
