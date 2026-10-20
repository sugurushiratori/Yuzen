# Numerical methods in Yuzen

This is the document to read when a calculation is slow, will not converge, or stops part of the
way through, and it is also the document to read before modifying the solver. It covers
discretization, the convection schemes, time integration, convergence tests, linear solvers,
parallel running, reproducibility, and every entry of the `Numerics` block. It is a shorter
companion to the user manual (`docs/manual/`), whose chapters "Discretization", "Time
integration", "Convergence and the correction floor", "Linear solvers", "Running in parallel",
"Bit-level reproducibility" and "Numerics and restart" give the full account; where the two
differ, the manual is authoritative.

* For the physical model and the meaning of the input parameters, see [`docs/theory.md`](theory.md).
* For verification results (comparison with analytical solutions, order of convergence, conservation error), see [`docs/validation.md`](validation.md).

---

## 1. Discretization

### 1.1 The staggered grid and the "site types"

A uniform rectangular grid is used. Given the grid lines `x = xmin + ix*dx`, there are only two
places along each axis at which a quantity can be stored.

```
0 : on the grid line               (x = xmin + ix*dx)
1 : halfway between two grid lines (x = xmin + ix*dx + dx/2)
```

The list of site types is therefore exactly three bits long and has exactly eight entries. Those
three bits are `Stagger<SX,SY,SZ>`, carried as a type (`src/grid/Stagger.h`).

| Bits (x,y,z) | Name | Position |
|---|---|---|
| 0,0,0 | `Node` | the grid point |
| 1,0,0 | `Xedge` | the edge along x |
| 0,1,0 | `Yedge` | the edge along y |
| 0,0,1 | `Zedge` | the edge along z |
| 0,1,1 | `Xface` | the face with normal x |
| 1,0,1 | `Yface` | the face with normal y |
| 1,1,0 | `Zface` | the face with normal z |
| 1,1,1 | `Cell` | the cell center |

**Moving half a cell along an axis flips that axis's bit.** Every change of site type is expressed
by that alone, and the accompanying movement of the index follows from a one-line rule: moving
+half a cell raises the index by one when the bit was `1`, and moving -half a cell lowers it by one
when the bit was `0`. There is no need to tabulate the offsets by hand.

In 1D only the x bit carries meaning, and in 2D the x and y bits do. A field `Field<GRID,S,T>`
takes the site type as a template parameter, so a stagger mismatch is a compile error. The
overloads of `C_FiniteDifferences::GradX` and its relatives, one per site type, exist for the same
reason. The unknowns sit on `Xedge` in 1D, on `Zface` in 2D, and the concentration of the `CZ`
model sits on `Cell`.

### 1.2 Ghost layers

Outside the boundary the grid carries **two ghost layers on each side** (`C_Grid*::NumGhostX()`).
This is because the thickness equation is fourth order, so evaluating it at a boundary point
reaches two neighbors away. The input `Grid.Nx` counts interior points only, so the total number
of points the grid holds is `Nx + 4`. `CheckGhost(s)` sorts every point into `Inside` (inside the
computational domain; this is the only place where the physical terms add their contributions),
`bXmin` / `bXmax` and the like (the boundary itself, where the boundary conditions fix the values;
in 2D the corners are added, and in 3D the edges as well), and `Outside` (the outermost ghosts).

### 1.3 From a point to a row number

The number of rows is fixed by the **total number of points including the ghosts**,
`imax() - imin() + 1`. The map from a point to its row number is collected in one place per
dimension (`Index()` in `GridTraits.h`).

| Dimension | Row number |
|---|---|
| 1D | $`\mathrm{row} = i_x`$ |
| 2D | $`\mathrm{row} = n_x\,i_y + i_x`$ |
| 3D | $`\mathrm{row} = n_x n_y\,i_z + n_x\,i_y + i_x`$ |

Ghost rows are marked as boundary-condition rows with `CRS_Matrix::SetAsBC()`. A marked row is
excluded from the time integration (`MakeMatrix_IA` / `MakeRHS` skip it), and the physical terms
consult `IsBC()` and do not add to it either. The coefficients of a boundary-condition row are
written once at initialization and survive every subsequent zeroing of the matrix.

The stencil is five points in 1D (`M m o p P`, reaching two neighbors away on either side of the
point itself) and 13 points in 2D (the 9-point cross plus the four diagonals `mm, pm, mp, pp`; the
diagonals are needed because the biharmonic operator contains a cross term in
$`1/(\Delta x^{2}\,\Delta y^{2})`$). The order in which the columns are inserted into the matrix is
fixed by the ordering of `T::FullIdx(k)` and its relatives (§8).

### 1.4 The mobility at a face

The film-height equation is a fourth-order *degenerate* diffusion: the mobility $`h^{3}`$ vanishes
where the film does, and the way a face mobility is built from the two cells either side decides
whether the scheme can go negative. The well-mixed model uses the positivity-preserving form of
Diez and Kondic (2002),

```math
D_{k}^{(x)} = \frac{2\,h_{k}^{2}\,h_{k+1}^{2}}{h_{k} + h_{k+1}},
```

the general prescription $`D_i = (h_{i+1} - h_i)/(g_{i+1} - g_i)`$, $`g = \int dh/D(h)`$, at
$`D = h^{3}`$ (`NonlinearTerm()` of `LaplacePressureKernel.h`). It vanishes as $`O(h^{2})`$ when one
side thins and is exactly zero when one side is dry; the arithmetic mean $`(h_k^3 + h_{k+1}^3)/2`$
would draw liquid out of a dry cell and drive $`h`$ negative. **Do not simplify it**;
`tests/test_positivity_mobility.cpp` pins it down.

**One face mobility for the three terms of the potential.** The disjoining pressure and the
vertical gravity take the same face mobility (`Mobility::Face`, the bits of the capillary one; 0
when either side is dry or is the ghost $`-h`$ of a `Pinning` edge), so that the flux of a face is
one mobility times the difference of the whole potential,
$`Q_f = D_f\,\Deltaarphi/\Delta x`$, $`arphi = \hat\sigma
abla^2(H+S) + \mathrm{Dp}\,\hat\Pi - \mathrm{Bo}(H+S)`$
(for $`\hat\sigma`$ uniform). Then the dissipation of the scheme,
$`\sum_f D_f(\Deltaarphi)^2/\Delta x`$, is a sum of squares, the discrete free energy never rises,
and a steady film is one with $`arphi`$ uniform to rounding: the discrete equilibrium is the minimum
of the discrete energy. The vertical gravity's term over a relief, $`D_f\,\Delta s`$, goes into
$`A`$ as $`K\,(h_a + h_b)\,\Delta s`$, $`K = D_f/(h_a + h_b)`$, the same share for the two cells.
Found through `Output.Energy` and corrected: the two took the mean $`(h_a^3 + h_b^3)/2`$, with which
the dissipation was not sign-definite: on the film of `spinodal_1d_unstable` grown to a steady film
it went negative and the free energy rose ($`5.9\times10^{-5}`$ in the measure of `analytic`
`energy_saturate_mono`, which now reads 0), and a mobility of $`h^3/2`$ survived against a dry
neighbor.

**The discrete energy** (`[Output.Energy]`, `model/well_mixed/EnergyBudget.h`). The free
energy is summed so that its derivative with respect to a cell's $`H`$ is exactly minus that cell's
potential (the gradient term over the faces, $`\tfrac12\hat\sigma_f(\Delta f/\Delta x)^2`$, the
others over the cells), and the dissipation is the power of the potential against the scheme's own
face fluxes, $`\mathcal D = \sum_i \varphi_i R_i\,\Delta A`$ (`HeightEquationBase::FluxPower`).
Then the semi-discrete system has $`d\mathcal F/dt = -\mathcal D`$ exactly, and Crank–Nicolson
$`(\mathcal F^{n+1}-\mathcal F^{n})/\Delta t + (\mathcal D^{n}+\mathcal D^{n+1})/2 = O(\Delta t^{2})`$:
the `Balance` column, the step's own (so that a run stopped and continued gives every row the same
value). Summed over a run, $`\mathcal F(T) - \mathcal F(0) + \int\mathcal D\,dt`$ falls at order
2.00 in Δt on leveling (periodic, with slip, between Symmetry walls with vertical gravity, 2-D) and
1.97 through a dewetting rupture (`analytic`, `energy_*`), and $`\mathcal F`$ never rises.

**With Navier slip** (`[Physics.Slip]`, theory §1.1) every mobility is taken from
`yuzen::Mobility` (`model/Mobility.h`); without it the capillary pressure, the body forces and the
Marangoni flux run their own no-slip expressions, and a run without the block is bit-identical to
one without slip in the code at all. With $`D(h) = h^{3} + a h^{2}`$, $`a = 3\Lambda`$
(the $`1/3`$ stays in the coefficient of the term),

```math
g(h) = \int\frac{dh}{D(h)} = -\frac{1}{a h} + \frac{1}{a^{2}}\ln\frac{h+a}{h}
     = -\frac{x - \ln(1+x)}{a^{2}},\qquad x = a/h,
```

and the capillary face mobility is $`D_f = \Delta h/\Delta g`$ again. Two cancellations are
avoided, and the switch points and their errors are:

* $`x - \ln(1+x)`$ for small $`x`$ (a slip length short against the film): below $`x = 0.1`$ it is
  summed as its series $`\sum_{n=2}^{18} (-1)^n x^n/n`$, whose remainder is below
  $`0.1^{17}/19 \approx 5\times10^{-19}`$ of the leading term; above, `x - log1p(x)` loses about
  $`2\epsilon/x \le 4\times10^{-15}`$ relative to the cancellation.
* $`\Delta g`$ for two nearly equal heights: where $`|h_2 - h_1| < 10^{-3}\,\bar h`$
  (`Mobility::kSeriesRel`), $`D_f = D(\bar h)/(1 + \delta^{2}K/24)`$ with
  $`\delta = h_2 - h_1`$ and $`K = D\,(1/D)''`$, whose truncation error is
  $`O((\delta/\bar h)^{4}) < 10^{-12}`$; above it $`\Delta h/\Delta g`$ loses about
  $`\epsilon\,\bar h/\delta < 10^{-12}`$ to the cancellation.

Against a reference without cancellation (Gauss–Legendre on $`\int dh/D`$ in long double) the
largest relative error is $`4.1\times10^{-13}`$ over $`\Lambda = 10^{-6}`$ to 3, heights
$`10^{-3}`$ to 7 and separations 0 to 100 $`h`$, both switches crossed; across each switch the
value is continuous to $`10^{-12}`$ (`test_positivity_mobility.cpp`). $`D_f = 0`$ on a dry side
still, but near it $`D \sim a h^{2}`$, so $`D_f`$ vanishes as $`O(h)`$, not $`O(h^{2})`$: the
degeneracy that keeps the film positive is weaker with slip. The disjoining pressure and the
vertical gravity take the same $`D_f`$; $`D`$ itself is extended to $`h < 0`$ as the odd
function $`h\vert h\vert(\vert h\vert + a)`$, as $`h^{3}`$ is. The upwinded body forces take
$`D(h_f)`$ of the face height, the Marangoni flux $`(h_f + 2\Lambda)\,h`$.

The `CZ` model advances the film height in the same conservative form,
$`\partial h/\partial t = -\partial Q/\partial x - E`$, with the same positivity-preserving $`h`$
factors; the flux of a face is the depth integral of the resolved velocity profile, so the mobility
carries the mean of $`1/\mu(c(\zeta))`$ weighted by $`3(1-\zeta)^2`$, the Marangoni flux the mean
weighted by $`2(1-\zeta)`$, and the surface tension is that of the surface concentration
(`FaceMobility()`, `FaceInvViscosity()`, `FaceSurfConc()`). A uniform viscosity gives
$`h^3/3\mu`$ exactly, which is why the two models agree where the concentration cannot affect the
film (`concz_vs_wm`: peak thickness to $`10^{-6}`$ relative, volume conserved exactly). Compare
them at a fixed step and with an output interval well below `EndTime`: adaptive steps put their
output times at different instants.

### 1.5 The convection schemes

Three kinds of term are convective: the advection of the concentration, and, in the well-mixed
film-height equation, the in-plane body forces (in-plane gravity and the centrifugal term) and the
Marangoni flux, which enter as a film height carried by a velocity rather than as the gradient of
a mobility. Each takes the value of a field at a face from the four cells about it. The scheme is
chosen by `Numerics.Convection.Concentration` (default `"UTOPIA"`) and
`Numerics.Convection.Height` (default `"Koren"`; refused with `CZ`, whose film-height equation
upwinds nothing); the keys are in §9.5.

| Scheme | Face value | Notes |
|---|---|---|
| UTOPIA | third-order upwind, $`c_1 = (-1,7,7,-1)/12`$, $`c_2 = (-1,3,-3,1)/12`$ | leading error $`-\vert u\vert\Delta x^3/12\,\partial^4 f`$; most accurate on smooth profiles |
| KK (Kawamura–Kuwahara) | UTOPIA with $`c_2`$ three times larger | three times the dissipation |
| Koren, VanLeer | limited: $`f_U + \tfrac12\Phi(a,b)`$ from the upwind side | bounded between $`f_U`$ and $`f_D`$; a smooth extremum is clipped to first order locally |
| Central4 | no upwind part; face value fourth order, divergence second | dissipates nothing |
| QUICK, Upwind1, Upwind2, Central2 | the classical forms | |

Sources: UTOPIA — B. P. Leonard, M. K. MacVean and A. P. Lock (1993), *Positivity-preserving
numerical schemes for multidimensional advection*, NASA TM-106055 (ICOMP-93-05); KK —
T. Kawamura, H. Takami and K. Kuwahara (1986), Fluid Dyn. Res. **1**, 145–162; Koren —
B. Koren (1993), in *Numerical Methods for Advection–Diffusion Problems*, Notes on Numerical Fluid
Mechanics **45**, Vieweg, 117–138; van Leer — B. van Leer (1974), J. Comput. Phys. **14**,
361–370; QUICK — B. P. Leonard (1979), Comput. Methods Appl. Mech. Eng. **19**, 59–98. Full
entries are in the manual's Further reading.

Koren's limiter is UTOPIA wherever $`0.4 \le a/b \le 4`$ and bounded elsewhere; nothing is divided
by a gradient that may vanish, and the mirror image of a stencil gives the same value to the last
bit. In the concentration equation a limited face value is frozen at the latest iterate so that the
frozen matrix is an M-matrix; with Crank–Nicolson **the bound holds for a Courant number up to one
only**. In the film-height equation $`A`$ holds the face weights frozen at the current film and the
Jacobian adds the rest of the derivative of the limited face height, so Newton keeps its quadratic
convergence (one-sided at the kinks of Koren's $`\Phi`$). The Marangoni flux is upwinded by the sign
of the flux itself, $`\sigma_c\,\Delta c`$, not of $`\Delta c`$ alone.

**Wall faces.** At the face of a `Symmetry` or `Pinning` edge the three upwinded terms of the
film-height equation carry no flux: their mobility there is zero in $`A`$, in the Jacobian and in
the velocity (`UpwindKernelBase::WallFace`). Their flux is not a gradient, and a mirrored ghost
does not make the face value zero.

**The Courant check.** With `"Koren"` or `"VanLeer"` in use the solver computes, after every step,
the Courant number of the step — $`\Delta t`$ times the larger of the rates of the states at its
start and end, a rate being the largest over the cells of
$`\sum_\text{axes}\max(\text{speed on the two faces})/\Delta`$ — from $`\vert u\vert`$ of the
advecting velocity (well-mixed concentration), from the layer volume fluxes $`\vert Q\vert/(h\Delta x)`$
and $`\vert W\vert/(h\Delta\zeta)`$ (`CZ` concentration), and from the characteristic speed
$`\partial Q/\partial h`$ of the upwinded fluxes, $`3\vert u\vert`$ of the body forces plus
$`2\vert u\vert`$ of the Marangoni flux (film height). The first step above one is reported in
`Warnings.txt` (grep for `Courant number`), with the count of such steps and the largest value at
the end of the run; the largest value is written to the log and to `JobInfo.txt` in any case. The
step is neither shortened nor rejected; to stay below one, shorten it (`Numerics.InitialDt`, or the
cap `MaxDtFactor`). Of the shipped cases, `nozzle_1d_rotating` (4.6) passes one in the film-height
equation, and `stress_all` reaches 0.94.

**0 ≤ c ≤ 1 without cutting or making solute.** With the evaporation limited to the solvent there
is, only rounding and the over- and undershoot of an unbounded scheme at a sharp front (KK,
UTOPIA, Central4, QUICK; Koren and van Leer above Courant 1) leave [0, 1]. In the well-mixed model
what a real cell holds above 1, $`h(c-1)`$, goes to its face neighbors in proportion to their room
$`h(1-c)`$, in up to four passes that are independent of the order of the cells, the rest spread
over every cell with room; a cell below 0 is filled the same way the other way round, what it lacks,
$`-hc`$, taken from its face neighbors in proportion to the solute $`hc`$ they hold, then from every
cell holding solute. The ghosts are formed again from the moved cells
(`HeightEquationBase::RedistributeExcess`). The sum of $`hc`$ is kept to rounding, and where no cell
leaves [0, 1] nothing is touched. Only a film holding more solute than liquid in all, or less than
none, would have a remainder cut or made up, reported in `Warnings.txt` (grep for `cut at c = 1`,
`made up at c = 0`). `JobInfo.txt` gives the solute moved over the run (`Solute moved from c > 1
[-]`, `Solute moved to c < 0 [-]`) and the solute cut and made up (`Solute cut at c = 1 [-]`,
`Solute made up at c = 0 [-]`). The working field the coefficients are taken from is capped at 1
within the Newton passes, but the limiter of Koren and van Leer is frozen at the iterate as solved:
frozen at the capped values, the rows on either side of a face would take different face values and
the advection would lose its flux form. `CZ` moves the excess of a cell to the cell below it in its
column; a cell below 0 is filled from the cells of its column holding solute, and a column whose
mean is below 0 from the columns next to it, in proportion to their solute (each scaled as a whole),
then from every column holding solute (`CZ_HeightEquation_2D::RedistributeDeficit`; it reports in
`JobInfo.txt` as the well-mixed model does). Left below 0, the undershoot of a front from 0.7 to
0.005 passed the divergence test of −0.1 and stopped the run (well-mixed QUICK and Central4, UTOPIA
at Courant 1.5; `CZ` UTOPIA at 0.4; found, corrected).

The defaults were chosen by the comparisons of the manual's section "Choosing the convection
scheme" (summarized in [`validation.md`](validation.md) §11): Koren for the film height because it
keeps the face height between its neighbors and survives steep fronts over a thin precursor, UTOPIA
for the concentration because it is the most accurate on smooth profiles. For a concentration with
sharp fronts Koren keeps $`c`$ bounded at the cost of accuracy at smooth extrema. Central2 and
Central4 dissipate nothing: where the model itself is unstable at the shortest resolved waves — in
`CZ`, a film stratified with $`\sigma_c\,\partial c/\partial\zeta > 0`$ and a weak diffusivity —
they let those waves grow.

### 1.6 The `CZ` concentration in conservative form

The `CZ` concentration equation is a finite-volume balance of $`hc`$ on the $`\zeta`$ grid. The
film-height equation stores the flux through each of its faces; the concentration equation splits
each into layer fluxes (the cell mean of the resolved velocity profile, adjusted so that the layers
of a column add up to the film-height flux face by face), and the flux through the $`\zeta`$ faces
follows from the discrete conservation of the grid, built up from the substrate, with the
evaporating solvent leaving through the surface face. The ghosts above the surface and below the
substrate make the total flux there — advection and the whole diffusive flux, its part along the
slope included — zero. The step is solved for the change of $`c`$, so a uniform $`c`$ gives exactly
no change, and the small solute residual of the linear solver is put back as a uniform shift of the
change, so that the solute of the step balances to rounding. The scheme of
`Numerics.Convection.Concentration` is used on every face, $`x`$, $`y`$ and $`\zeta`$ alike; the
limited schemes are re-solved until the limiter has settled at the midpoint of the step.

### 1.7 The well-mixed concentration in conservative form

The well-mixed concentration equation is the finite-volume balance of the solute $`hc`$ of each cell
over the step (`ConcEquationBase.h`, `ConcStep.h`):

```math
h^1 c^1 - h^0 c^0 = -\Delta t \sum_f s_f \left[\theta G_f^1 + (1-\theta) G_f^0\right]
  + \Delta t\, S\, c_S - \Delta t \left[\theta r^1 c^1 + (1-\theta) r^0 c^0\right],
\qquad G_f = F_f c_f - D(c_f)\, h_f \frac{c_\mathrm{hi} - c_\mathrm{lo}}{\Delta x},
```

$`s_f = \pm 1/\Delta x`$, $`\theta = 1/2`$. $`F_f`$ is the volume flux the film-height equation puts
through the face, summed over the enabled terms from what each holds in $`A`$ for it (`FaceFlux` of
the kernels: capillary, vertical gravity, disjoining, and the upwinded in-plane body forces,
centrifugal and Marangoni fluxes with their face heights); $`c_f`$ the face value of
`Numerics.Convection.Concentration` upwinded by the sign of $`F_f`$; $`S`$ the nozzle's supply, at
the middle of the step as the film-height equation takes it, and $`c_S`$ the concentration it arrives
at: $`\theta c^1 + (1-\theta) c^0`$, the concentration it lands on, with
`Physics.Source.ConcType = "local"` (the default), or `Physics.Source.Conc` with `"supplied"`, the new
liquid mixed through the column (with no flow, $`c^1 = (h^0 c^0 + \Delta t S c_S)/(h^0 + \Delta t S)`$
exactly); $`r`$ the 1-D sideways outflow. Evaporation takes
solvent only and is in $`h^1`$. The unknown is the change $`d = c^1 - c^0`$: subtracting $`c^0`$
times the film's own step leaves an equation whose right-hand side is formed from face values
$`c^0_f - c^0 = \sum_k w_k (c^0_k - c^0)`$ and differences, exactly zero for a uniform $`c^0`$ without
evaporation (and $`\Delta t\,E\,c^0`$ with it, which gives $`h^1c^1 = h^0c^0`$ where nothing else acts)
and without liquid supplied at another concentration (which adds $`\Delta t\,S\,(c_S - c^0)`$; with
`"local"` the supply adds $`\theta\,\Delta t\,S`$ to the diagonal instead). `History.dat` counts the
nozzle's solute (`SuppliedResin`) from the same supply term: with `"supplied"` it is
$`\sum \Delta t\,S\,c_S`$ over the cells and steps, so `ResinVolErr` of a closed film with a nozzle is at
rounding in both modes (`nozzle_solute` 1.3e-13, `nozzle_supplied` 5e-14).

On every Newton pass the equation is solved with the film of the pass and its fluxes. Once the
iteration has converged the step is completed (`AfterNKLoop`): the film is the one it converged
to; the fluxes, outflow and supply are taken there (with the concentration of the last pass, and the
evaporation as that pass took it); what their divergence misses of the step the film took, cell by
cell — the residual of the step there, less its mean, which is the volume error and stays as that —
is put into the inner faces as running sums (along $`y`$ each row's sum spread over the row, then
along $`x`$ within each row; `ConsistentFluxes`), so that the face fluxes are exactly those of the
step; no edge face and no wall is touched. (The film the fluxes give, $`h^0`$ less their divergence,
is not kept instead: on a film far thicker than its cells the residual is at the rounding of
enormous fluxes, §4.1, and that film goes negative in `nozzle_1d_rotating`.) The concentration is
solved once more, and the solute of each cell formed from the balance with each face's solute flux
formed once (`Reconstruct`): $`\sum hc`$ changes by the edge fluxes and the sources only, to rounding,
whatever the tolerance of the linear solve. Those edge fluxes are what `History.dat` counts as crossing
an open edge (`FreshInlet`, `Outlet`): the volume flux of each edge face of an open side, at the start
of the step and at the film it ended on, weighted as the step weights them (`EdgeFluxRate`; an edge
face is left alone by `ConsistentFluxes`), and the solute flux `Reconstruct` formed there, carried and
diffused (`EdgeOutflow`), summed over the committed steps only (`EdgeOutflowVolume`), the totals kept
in the state file. `VolErr` and `ResinVolErr` of an open film are then at rounding as on a closed one
(`inlet_outlet`, `solute_balance`: $`10^{-15}`$ to $`10^{-14}`$). Counted from the velocities and the
face values of the field after each step, integrated by the trapezoidal rule, they were found off by
the error of the scheme, $`6\times10^{-7}`$ on 100 cells, falling at order 1.8 in $`\Delta x`$;
corrected. A limited scheme is re-solved with its limiter frozen at
its own answer until that settles ($`10^{-12}`$, at most 8 times). Then $`c`$ is brought into
[0, 1] (§1.5). A pass, or the end of a step, whose film has a cell with $`h \le 0`$ (a step far too
long for the film) leaves $`c`$ as the pass before had it, and the divergence test stops the run on
the film. That test also takes the lowest $`c`$ of the start of the step and of the solve, before a
cell below 0 is filled, so a fill cannot hide a solve gone far below 0 or a field file read with
one.

Solved for $`c`$ in non-conservative form, a front from 0.9 to 0.3 on a closed periodic film
carried by the flow and drying by Raoult's law lost $`1\times10^{-5}`$ to $`3\times10^{-4}`$ of its
solute over 0.3 of a crossing (50 to 200 cells, six schemes); found, corrected: $`4\times10^{-16}`$
to $`3\times10^{-15}`$ in 1-D, $`2\times10^{-14}`$ in 2-D (`tools/analytic.py conc_drying`). A uniform
$`c`$ under a strong flow, a nozzle, gravity, capillarity and (1-D) spin stays uniform to
$`6\times10^{-16}`$ (`conc_uniform`).

## 2. Time integration

### 2.1 The form that is solved

Every equation is rearranged into the following form before it is solved
(`numerics/ImplicitSolver.h`).

```math
\frac{\partial\varphi}{\partial t} = -A\,\varphi - B
```

$`A`$ is the part proportional to the unknown (a matrix) and $`B`$ is everything else (a vector).
Time is advanced with the $`\theta`$ method, where $`\theta = 0 / 0.5 / 1`$ correspond to forward
Euler / Crank-Nicolson / backward Euler. Both the thickness equation and the concentration
equation are built with **Crank-Nicolson ($`\theta = 0.5`$)**.

```math
(I + \theta\,\Delta t\,J)\,\delta\varphi = \varphi_\mathrm{old} - \varphi_\mathrm{guess} - \Delta t\left[ \theta\,A\cdot\varphi_\mathrm{guess} + (1-\theta)\,A\cdot\varphi_\mathrm{old} + B \right]
```

> **★ A term in $`B`$ has to be averaged in time as well ★**
> The $`\theta`$ weighting applies to $`A`$; a term in $`B`$ gets it only if $`B`$ is formed as the
> same average. The well-mixed thickness equation does so for its two terms in $`B`$: the
> evaporation as $`(1-\theta)B(\varphi_\mathrm{old}) + \theta B(\varphi_\mathrm{guess})`$, the liquid supply at the
> middle of the step. Taken at the start of the step alone, $`B`$ is first order in dt. A new
> physical term goes into $`A`$ if it can be read as a coefficient times the unknown -- the
> disjoining pressure $`\hat{\Pi} = (\hat{\Pi}/H)\cdot H`$ is the example -- or gets the same average in $`B`$.

**How the product with the matrix is formed.** For the well-mixed thickness equation the
products $`A\cdot\varphi`$ are not all formed from the entries of $`A`$ (`ApplyA` in
`HeightEquationBase.h`). The capillary pressure and the vertical gravity form theirs **face by
face** (`AddConservative` of their kernels): the flux through a face is its mobility times a
difference of the film taken between neighboring values, and the two cells that share the face
take it from the same values in the same order, so one loses exactly what the other gains. Formed
from the entries, which are of order $`D/\Delta x^{4}`$, the sum over the cells would carry
$`\varepsilon_\mathrm{mach}\Vert A\Vert\vert\varphi\vert`$ per cell — a volume error of
$`2\times10^{-7}`$ of the film over `nozzle_1d_rotating`. The in-plane body forces, the Marangoni
stress, the disjoining pressure, the sideways outflow of the 1-D centrifugal term and the rows of
the boundary conditions keep the product with their entries. **The boundary rows are satisfied
exactly** after every linear solve (`CRS_Matrix::SatisfyBCRows`), each set of rows that refer to
one another together, so a `Symmetry` ghost is the exact mirror of the layer inside, a `Periodic`
one the exact copy of the other end, and no flux leaks through an edge face.

### 2.2 The Newton-Kantorovich outer iteration

The thickness equation is non-linear because its mobility $`h^{3}/(3\mu(c))`$ depends on the
unknown, so Yuzen iterates with Newton's method within a single time step until the non-linearity
itself is solved (`numerics/NewtonKantorovich.cpp`). One iteration has five stages.

```
1. EveryNKLoop()     update the fields from the guess φ_guess and solve the coupled side (concentration)
2. SetMatrixA()      assemble the matrix A
3. SetMatrixJacobi() assemble the Jacobian (once every JacobianUpdateInterval iterations)
   MakeMatrix_IA()   form I + θ dt J
4. linear solve      obtain the correction δφ
5. φ_guess += α δφ   (α = 0.5 only when the correction has grown since the previous iteration)
```

The iteration starts from the film extrapolated from the last two steps,
φⁿ + (dt/dt_old)(φⁿ − φⁿ⁻¹), unless that would leave a cell that has liquid without any; an attempt
from it that fails is made once more from φⁿ at the same step, before the step is cut or the run
stopped (a nozzle that turns on is not in the last two steps). Each pass's linear solve is solved
only as far as the iteration needs (§7.1, Eisenstat-Walker). The two cut the Newton iterations by up
to a third and the linear iterations by up to 38 % on the cases measured; the convergence tests
stay relative to the first correction of the step (§4.3), so the converged answer stays where it was.

The concentration equation, and the thickness equation of the `CZ` model, are weakly non-linear,
so they use the semi-implicit path (`numerics/SemiImplicit.cpp`), which freezes the coefficients and
solves a linear system, with no outer iteration. Frozen at the start of the step that is first order
in dt, so `CZ` solves each step twice -- a predictor with the coefficients of the start, then the
same system from the same start with every coefficient and every term of B at the midpoint (the
implicit midpoint rule, second order), the thickness equation with the concentration extrapolated
to the midpoint from the last two steps. The velocity the solute is carried with is formed from the
film halfway to the end state. Without the Marangoni term the concentration equation is then solved
once; with it, the film and the concentration are solved again at the midpoint concentration the
last pass found, with Aitken relaxation, until that midpoint moves by less than a tenth of the
extrapolation's miss (at most eight passes; `JobInfo.txt` reports how many were taken). On a flat
film the gain of a plain pass exceeds 1 only where a wave of the model grows faster than 2/dt (the
gain is derived and measured in [`validation.md`](validation.md) §7.3); a step whose passes do not settle
within eight has failed, as one whose linear solve does not converge: the time-step control cuts it,
and at a fixed step the run stops. `YUZEN_DEBUG_COUPLING` prints every pass. The well-mixed concentration equation is
solved inside the Newton loop on the current guess of the end state, so its Crank–Nicolson average
is already that of the two ends.

### 2.3 The Jacobian

Each physical term adds its own analytic Jacobian through `SetJacobiMatrix()` (the contract in
`PhysicsTermBase.h`). The difference coefficients are built once by `SetLinearCoef()` from the mesh
and the scheme alone, while the field-dependent part is separated into `NonlinearTerm()` and its
derivative `NonlinearJacobi()`.

**The Jacobian is allowed to be approximate.** It changes only how fast the iteration converges,
not where it converges to (the fixed point of the iteration). What fixes the fixed point is the
residual `Rhs`, not the Jacobian. Within the thickness equation this property is exploited in one
place only, the modified Newton method below; the terms themselves are exact, and
`tests/test_jacobian.cpp` checks each against a finite difference, in 1-D and 2-D and with every
convection scheme.

* The Jacobian starts as a copy of A. For the gravity and centrifugal terms the gradient side is
  linear in h, so A already carries its derivative and the term adds only that of the mobility.
* `Numerics.HeightEquation.JacobianUpdateInterval` greater than 1 turns the scheme into a modified
  Newton method, in which the Jacobian is rebuilt only once every few iterations (`_crsIA` can be
  held with it, since it depends on nothing but `_crsJacobi` and `dt`). The 1D breakdown of one
  iteration is `SetMatrixA` 19 % / `SetMatrixJacobi` 43 % / `IA + RHS` 2 % / the linear solve 37 %,
  so the Jacobian assembly is the single largest item.

**A gradient side that is not linear needs its own derivative.** The disjoining pressure goes like
$`\Pi \sim h^{-n}`$, is extremely stiff, and is written into A as $`(\hat\Pi/h)\,h`$, so its term
adds the mobility derivative and $`\hat\Pi' - \hat\Pi/h`$ on the gradient side. Adding all of
$`\hat\Pi'`$ would count $`\hat\Pi/h`$ twice and leave the Jacobian 50 % off; with the correct
derivative the disjoining-pressure cases take 30–55 % fewer Newton iterations.

## 3. Coupling of thickness and concentration

Thickness and concentration are solved as a block Gauss-Seidel iteration inside the outer Newton
loop.

```
for each time step, k = 0,1,2,...
  form the face volume fluxes F of h^k
  solve the concentration equation over dt  (F, coefficients and end film at h^k)  → c^{k+1}
  cap the c of the coefficients at 1        (ClampConc)
  take one Newton step on the thickness     (Jacobian is ∂F_h/∂h only)            → h^{k+1}
once converged, to h^1:
  the fluxes at h^1, made those of the step h^0 → h^1 (ConsistentFluxes)
  solve the concentration equation once more                                     → c^1
  bring c^1 into [0, 1] without cutting solute (RedistributeExcess)
```

The off-diagonal blocks ($`\partial F_h/\partial c`$, $`\partial F_c/\partial h`$) are dropped for
exactly the reason given in §2.3: they slow the convergence down but **do not move the converged
answer**. If the loop converges to $`\delta h \to 0`$ together with $`c^{k+1} = c^{k}`$, then the
thickness equation and the concentration equation are satisfied simultaneously, and that state is
precisely the fully implicit coupled discretization of the two. For that reason **the convergence
test is made on the larger of the two corrections**,
$`\max( \vert\delta h\vert_\infty , \vert\Delta c\vert_\infty )`$.

> **★ Do not test on the thickness correction alone ★**
> $`\vert\delta h\vert`$ can become small while the concentration is still moving, so a test that
> leaves the concentration out lets the iteration exit early. What it exits to is a state that
> satisfies the thickness equation only, not the fully implicit coupled solution. Propping the
> iteration up with a minimum iteration count is to be avoided, since it decides the accuracy with
> a number that has nothing to do with the problem.

The concentration correction is measured **after the clamp (`ClampConc`)**, because the clamped
value is the one that enters the coefficients of the physical terms and therefore the one that
decides whether the coupling has settled. A model that does not report a coupling (`CZ` solves the
thickness equation semi-implicitly too and has no outer iteration at all) leaves the reported value
at zero, so there the test reduces to the thickness. `NewtonMinIterations` can force a minimum iteration
count (default 1), but since the concentration correction is part of the test the default of 1 is
normally the right value. To see whether it is the thickness or the coupling that is holding
convergence up, set the environment variable `YUZEN_DEBUG_RESIDUAL`; each outer iteration then
prints a line to standard error (`[resid] iteration number, thickness correction, concentration
correction`), and where the iteration stagnates (§4.3) a line `[stall] iteration, correction, first
correction, the floor 1000 eps |phi| it was judged against, 1 accepted / 0 failed`.

## 4. The convergence test and the "correction floor"

This is the most practically important section.

### 4.1 The floor of the correction

Newton's method converges quadratically, so the correction falls cleanly. Then it stops at a
certain height, and from there on it merely scatters instead of decreasing. On the spreading film of
the table below at `Nx = 2000`:

```
iteration 0 : 6.231e-02      <- quadratic convergence, falling cleanly
iteration 1 : 1.899e-03
iteration 2 : 1.494e-05
iteration 3 : 1.037e-09
iteration 4 : 1.931e-16      <- this is the "floor"
iteration 5 : 2.220e-16
iteration 6 : 2.189e-16      <- from here it wanders between 1.5e-16 and 2.2e-16
```

The floor is **cancellation in the assembly of the right-hand side**, and it has nothing to do with
the accuracy of the linear solver. The right-hand side is

```math
\mathrm{Rhs} = (\varphi_\mathrm{old} - \varphi_\mathrm{guess}) - \Delta t\left[ \theta\,A\cdot\varphi_\mathrm{guess} + (1-\theta)\,A\cdot\varphi_\mathrm{old} + B \right]
```

and $`A`$ is a fourth-order spatial difference, so its coefficients are of order
$`1/\Delta x^{4}`$: at `Nx = 2000` (a non-dimensional `dx = 1.72e-4`) the row sums reach `~6e15`.
Formed as the product of those entries with $`\varphi`$, $`A\cdot\varphi`$ is an answer of about
`1e5` obtained by canceling terms of that magnitude against one another, with an absolute error of
$`\varepsilon_\mathrm{mach} \times \Vert A\Vert \times \vert\varphi\vert \sim 1`$ multiplied by
$`\Delta t`$. Since the inverse of $`I + \theta\,\Delta t\,J`$ is a contraction, the floor of such a
product is roughly

```math
\varepsilon_\mathrm{mach} \times \frac{\Delta t}{\Delta t_\mathrm{explicit}}
```

```math
\Delta t_\mathrm{explicit} \sim \frac{3\mu\,\Delta x^{4}}{H^{3}}
```

proportional to how far beyond the explicit stability limit the step goes. The well-mixed
thickness equation does not form its stiff terms that way: the capillary pressure and the vertical
gravity give their $`A\cdot\varphi`$ face by face (`ApplyA` in `HeightEquationBase.h`,
`AddConservative` in the two kernels), the flux of a face its mobility times a difference of the
film, the differences taken first between neighboring values, and the two cells that share a face
take the same flux. What cancels is the variation of the film over the stencil rather than the
film, so on the films below the floor is that of $`\varphi`$ itself, a few
$`\varepsilon_\mathrm{mach}`$. The formula remains the bound for what is formed as a product (the
in-plane body forces, the Marangoni stress, the disjoining pressure, whose entries are of order
$`1/\Delta x`$ or $`1/\Delta x^{2}`$) and for a film whose variation over a stencil is comparable to
the film. `YUZEN_DEBUG_FLOOR=1` prints the floor.

| Case | $`\Delta x`$ [µm] | $`\Delta t`$ (non-dimensional) | $`\Delta t/\Delta t_\mathrm{explicit}`$ | bound $`\varepsilon_\mathrm{mach}\,\Delta t/\Delta t_\mathrm{explicit}`$ | measured floor |
|---|---:|---:|---:|---:|---:|
| spreading, Nx=500 | 1.2 | 3.21e-9 | 4.8e3 | 1.1e-12 | 1.1e-16 |
| spreading, Nx=2000 | 0.3 | 3.21e-9 | 1.2e6 | 2.7e-10 | 1.8e-16 |
| spreading, Nx=2000 | 0.3 | 3.21e-11 | 1.2e4 | 2.7e-12 | 1.1e-16 |
| flat film (`leveling_1d_reg`) | 3.33 | 3.21e-8 | 8.1e2 | 1.8e-13 | 5.5e-17 |

The spreading film is `stress_wetting` at the `Nx` given, run with the adaptive step for 1 ms (0.1 ms at
`Nx = 2000`) and then stepped at the fixed `dt` of the table; the floor is the median of the
corrections from the fifth iteration on.

* **No particular value is universal.** On these films the floor does not grow as
  $`\Delta x^{-4}\Delta t`$, but a term formed as a product, a film steep over a few cells or another
  case can put it higher, as high as the bound. Measure it rather than assume it.
* **Machine precision is the limit.** The correction does not fall below
  $`\varepsilon_\mathrm{mach}\,\vert\varphi\vert`$. The linear solver is not what sets it (the
  measured `bcg_resid ~ 1e-25`), and tightening its tolerance does not move the floor.

### 4.2 What an absolute-only test does

> **★ Do not use an absolute-only convergence test ★**

Where the floor is proportional to $`\Delta t`$, an absolute threshold below it creates the relation
"make $`\Delta t`$ smaller and the floor falls until the test is met". That relation closes into a
feedback loop (and a threshold below $`\varepsilon_\mathrm{mach}\,\vert\varphi\vert`$ is never met).

1. A `NewtonEpsilon` below the floor can never be satisfied.
2. The iteration has effectively converged after five iterations, but it spins on to `NewtonMaxIterations`.
3. Having hit the iteration cap, it reports a convergence failure.
4. The caller (`TimeAdvance`) sees the failure and **halves the time step**.
5. A smaller step lowers the floor, so the run locks onto the step at which `NewtonEpsilon` is just met
   and oscillates there (`dt-` and `dt+` appear alternately).

The finer the mesh, the higher a floor of that kind and the smaller the step the run locks onto: on a
mesh that resolves a contact line the loop can cost orders of magnitude in the physical time
reached.

### 4.3 The two tests that avoid the floor

**(1) The relative convergence test `NewtonEpsilonRel`**

The iteration is taken to have converged once the correction has fallen to `NewtonEpsilonRel` times
**the first correction of that time step** (default `1e-6`). The statement it makes is "the
non-linearity of this time step has been solved to six digits", which is independent of the scale
of the solution, of the time step and of the mesh. The floor and `dt` cancel between numerator and
denominator, which breaks the coupling of §4.2. On the history above the threshold is
`6.231e-02 × 1e-6 = 6.23e-08`, so the iteration exits at iteration 3. Setting it to `0` or below
switches the relative test off and leaves the absolute test alone. From the extrapolated start
(§2.2) the first correction is what the extrapolation missed, and the test stays relative to that:
relative to how far the extrapolation moved the film instead, the iteration stopped a pass earlier
and the volume of a nozzle film thickened a hundredfold drifted to |VolErr| = 5e-11
(`tools/analytic.py` nozzle_balance; 7e-14 so).

**(2) The stagnation test**

This looks at "how many iterations in a row fail to improve the best value so far". An iteration
that does not get below 0.9 times the smallest correction seen so far counts as one failure, and
more than 5 failures stops the loop. Measuring the change ratio
(`|resid - resid_old| / resid_old`) instead does not fire in a situation where the correction
scatters by tens of percent, or by several times, on the floor, which is why the best value is the thing to watch.

When the test fires, the iteration has converged as far as it can, and the step is accepted
silently, if the absolute test is met, or the correction has fallen four or more decades below the
first correction (the non-linearity of the step resolved to four digits), or it is within
`1000 eps |phi|` (`|phi|` the largest film value): the film cannot be updated more finely than
`eps |phi|`, and the floors of §4.1 are about that. Otherwise the iteration has
stopped **far from convergence** and the step has failed, like a linear solve that does not converge
or a run out of `NewtonMaxIterations`. With `Numerics.TimeStepControl.Enable = true` it is cut by
`ChangeFactor` and solved again, the first one reported in `Warnings.txt` and all counted in
`JobInfo.txt` (`Steps failed: the Newton correction stagnated`); a cut below `MinDt` aborts the run
as any other. At a fixed step the run stops and names the stagnation.

At the default tolerances every stagnation in the shipped cases and the checks is of the second
kind, and far into it: on the contact-line stress case (`stress_wetting`, the first step) the
correction stopped at 0.068 to 25 times the first, on the thick nozzle film with liquid supplied at its own
concentration (`nozzle_supplied`) at 8e-4 to 1.9 times; none fell between the relative test (1e-6)
and 8e-4. Kept, such a step is not a solution of the step: the nozzle film kept its residual as a
volume error of 2.2e-8 (found; corrected: cut, 2e-13), and the contact-line film at a fixed step
went negative and ended as diverged.

The four decades are deliberately not tied to `NewtonEpsilonRel`. With `NewtonEpsilonRel` at 1e-11 (`NewtonEpsilon` 1e-16)
the same nozzle film levels off at 1e-11 to 1e-4 of the first correction in two steps out of three
(the limited upwinding at Courant numbers up to 150 is not smooth), and those steps are accepted;
failing them would chase the relative test the way §4.2 chases an absolute one. The run cuts 33
steps that stagnated above 1e-4, takes 3358 steps and ends with `VolErr` 1.5e-9, against 2e-13 at
the default: what the accepted steps left. A tighter `NewtonEpsilonRel` buys nothing on a film whose
iteration cannot reach it.

### 4.4 What to look at when it does not converge

The input keys are `NewtonEpsilon` / `NewtonEpsilonRel` / `NewtonMaxIterations` / `NewtonMinIterations` in
`Numerics.HeightEquation` (§9.3).

1. Use `YUZEN_DEBUG_RESIDUAL` to see whether it is the thickness or the concentration that is
   holding the iteration up.
2. Use `YUZEN_DEBUG_FLOOR` to see where the correction stops falling.
3. If the floor sits above `NewtonEpsilon`, rely on the relative test rather than lowering `NewtonEpsilon`.
   Lowering `NewtonEpsilon` is the maneuver of §4.2 run in reverse.
4. If the warning "the Newton correction stagnated far from convergence" does appear (and
   `JobInfo.txt` counts the steps cut for it), then it genuinely is not converging and the floor is
   not the explanation. A non-smooth ingredient (upwind switching, or the clamping applied to
   evaporation and to the concentration) may be at work, or a step far too long for the state, as
   at the start of `stress_wetting`.

## 5. Time-step control

The step changes only when `Numerics.TimeStepControl.Enable = true`. It is decided by the thickness
equation; the concentration equation follows whatever the thickness equation asks for.

### 5.1 The local truncation error

The local truncation error is estimated from two consecutive time differences
(`C_ImplicitSolver::LocalTruncationCriterion`).

```math
\mathrm{Ek}_i = 2\,\frac{\Delta t}{\Delta t_0}\,\frac{\Delta t_0\,h^{+}_{i} + \Delta t\,h^{-}_{i} - (\Delta t_0+\Delta t)\,h^{0}_{i}}{(\Delta t_0+\Delta t)\,w_i}
```

```math
\mathrm{criterion} = 1 - \max\vert \mathrm{Ek}\vert / \mathrm{Criterion}
```

Here $`\Delta t_0`$ is the step last taken and $`\Delta t`$ the current one; for the first step,
$`\Delta t_0`$ is `InitialDt`, so the first step is judged too. A positive `criterion`
accepts the step; a negative one shrinks it. A step whose solve fails — the Newton iteration of the
well-mixed film (not converged in `NewtonMaxIterations`, or stagnated far from convergence, §4.3),
or a linear solve of the film or of the concentration in either model — is shrunk
in the same way and solved again; with the step control off it stops the run, naming the solve,
since a fixed step has no shorter one to be retried at. $`\Delta t_0`$ must be the step that produced
$`h^{-}`$: taking the step before its last change, or the rejected step on a retry, leaves a
first-order term in the estimate and holds the step down. `Criterion` is not the error itself: a
leveling film at `Criterion` $`= 10^{-5}`$ takes 69 steps for an error in its decay of
$`2.9\times10^{-4}`$, so choose it by the error you need. The weight $`w_i`$ that normalizes the
error is
$`\max( \vert h^{0}_{i}\vert , \mathrm{ScaleFloor} \times \max\vert h^{0}\vert )`$, so
`ScaleFloor = 0` (the default) gives a completely local relative error, dividing by the local
thickness; `1` gives something close to an absolute error, normalizing by the global maximum alone;
and values in between judge the body of the film relatively while refusing to over-weight the thin
regions.

Read as a sentence it says "anywhere thinner than `ScaleFloor` times the maximum thickness is not
judged more harshly than that", which also serves as the division-by-zero guard in a dry region
($`h \to 0`$). A precursor film ($`H_* = 0.02`$) is a fiftieth of the body of a drop
($`H \sim 1`$), so with the default a wobble inside the precursor film weighs fifty times as much
as the same wobble in the body. **This is still not the usual reason a contact-line calculation
sees its step collapse** (§4.2); suspect the convergence test first.

### 5.2 Growing and shrinking the step

| Direction | Condition | Operation |
|---|---|---|
| shrink | `criterion <= 0` | `dt *= ChangeFactor` |
| grow | `criterion > 0.9` holds for `KeepSteps` steps | `dt /= ChangeFactor` |

**A shrink acts at once, on the step being retried; a lengthened step applies from the next step,
and the clock advances by the step that was solved.**

The step is capped above and below. The upper bound is the smaller of
`MaxDtFactor` $`\times`$ (minimum mesh width)$`^{4}`$ and the smallest output interval; the lower
bound is `MinDt`. Both bounds apply with the step control off as well, so a fixed step is still
capped and floored; the `CZ` model applies neither, its step being capped by the output intervals
alone and floored at a built-in `1e-15`.

As soon as a request below the lower bound is made, the run aborts. A calculation that shrinks that
far does not have a time-step problem, so there is no point in lowering `MinDt` and carrying on.
The output interval caps the step because a step that straddles an output interval shifts the
output time. The input keys are the 7 entries of `Numerics.TimeStepControl` (§9.2).

**Field output at listed times.** `Output.FieldData.Times`, or a `LogSchedule` (`Start`, `End`,
`PerDecade`, read as Start·10^(i/PerDecade)), replaces the field interval, and then caps only the
step that would pass a listed time (`C_TimeSteps::ApplyListedTimes`): without `ExactTimes`, the
step is held to the gap between the two listed times around it and the field written at the step
closest to each; with `ExactTimes`, a step that would pass the next listed time is shortened to land
on it, and the step after it goes back to the step the run was on before. The shortened step is the
step taken (`DT().old`, `LastStepDT()`), so the error estimate of the next step differences over
the real distance; when the step control asks for less than it, or the step was retried shorter,
its request stands. `leveling_1d_logout` writes 141 times over seven decades in 403 steps (278
without `ExactTimes`), each within rounding of its listed time; held to the shortest gap it would
take at least 10^7. The other outputs still cap the step at their intervals.

## 6. Termination

At the end of every step `CheckTermination()` applies three tests in order. **(1) The specified
time or the specified step count has been reached**, in which case the run stops. The two are
reported separately (`Time reached to the EndTime.` for `Numerics.EndTime`,
`Step count reached Numerics.MaxSteps.` for `Numerics.MaxSteps`; `MaxSteps` defaults to 0, no
limit, so unless it is set the end time is what stops the run). **(2) The film has dried
out** (the concentration selected by `Physics.Evaporation.EndConcType` -- `"mean"`, `"max"` or
`"min"`; any other string is refused when the input is read -- has exceeded `EndConcValue`), in which
case the run prints `Dry out.` and stops. The third is the steady-state test.

> **★ Do not judge steadiness on the change per step ★**
> $`\max\vert h_\mathrm{new} - h_\mathrm{old}\vert / \max\vert h_\mathrm{new}\vert`$ becomes small
> simply because the time step is small. A calculation that resolves the contact line at
> `dt ~ 1e-11` would be declared "steady" within a dozen or so steps of the start and stop without
> complaint. The more finely the case is resolved the sooner this happens, which makes it hard to
> notice.

What the test actually uses is the relative change **per unit time**.

$`\mathrm{ChangeRate} = \max\vert h_\mathrm{new} - h_\mathrm{old}\vert / ( \max\vert h_\mathrm{new}\vert \times \Delta t )`$ `[1/s]`

$`\mathrm{ChangeRate} \times (\mathrm{EndTime} - \mathrm{StartTime}) < \mathrm{SteadyStateCriterion}`$
  = "how much the shape would change if it went on changing at the present rate for the whole span of the run"

The span runs from the start of the run (0, or `Restart.RestartTime`) to `EndTime`. This is a
dimensionless quantity that depends neither on the time step nor on the time scale `T0`.
The default `1e-3` means "at the present rate, the whole span of the run would move it by only
0.1 %". Once this has held for more than `SteadyStateSteps` consecutive steps the run prints
`Steady state.` and stops. The first
five steps are not judged at all (`ChangeRate` is reported as 1 for them). Setting
`SteadyStateCriterion <= 0` disables the test.

`ChangeRate` is written out as a column of `CalcInfo.dat`, so look at it first when choosing
`SteadyStateCriterion`. The input keys are in §9.1.

## 7. Linear solvers and parallel running

### 7.1 The backends

There are two implementations of the one job of solving a system of linear equations, and they
are interchangeable behind a single interface (`ILinearSolverBackend`). The calling code knows only
`LinearSolver`. The choice is made by `Numerics.HeightEquation.LinearSolver.Backend` and
`Numerics.ConcEquation.LinearSolver.Backend`, so the thickness equation and the concentration
equation can be given different backends.

| Value | Implementation | Character |
|---|---|---|
| `auto` (default) | amgcl for a 2-D film on a build with AMGCL, otherwise eigen | the one to pick when in doubt |
| `eigen` | Eigen `BiCGSTAB` + block-Jacobi `IncompleteLUT` | depends on Eigen alone; parallel. **Bit-identical whatever the thread count** |
| `amgcl` | AMGCL `bicgstab` + smoothed-aggregation AMG + SPAI0 | fastest in 2-D on most sizes. OpenMP-parallel including the preconditioner |

Asking for `amgcl` on a build without AMGCL emits a warning (`LinearSolver: Backend="amgcl" was
requested, but this build has no AMGCL. Falling back to "eigen".`) and uses `eigen`. Any other value
is refused (``Invalid value. Numerics.HeightEquation.LinearSolver.Backend must be one of [`auto',`eigen',`amgcl'].``, and the same for the
`ConcEquation`) and the run stops with exit code 1.
The `###### Solver Info ######` section of `JobInfo.txt` names the backend each equation actually
used, so it shows what `auto` resolved to and records a fall-back.

What each does (manual, sec:linsolve:what):

* **eigen, `ILU`: block-Jacobi ILUT.** The rows are split into contiguous blocks, one per 1024 rows
  and at most 128 (fewer than 2048 rows: one block, the ILUT of the whole matrix); each diagonal
  block is factorized by `IncompleteLUT` (fill 10, drop tolerance 1e-10), the blocks in parallel.
  The blocking follows from the size alone, so the result does not depend on the thread count. It
  takes 4 to 9 times the iterations of one ILUT of the whole matrix on the 2-D films measured and
  is 1.3 to 3 times faster on one thread all the same, 3.7 to 21 times at the best thread count.
* **amgcl, `ILU`: smoothed-aggregation AMG with SPAI0**, built anew at every solve. A system of more
  than 3000 unknowns is coarsened down to 500 (AMGCL's default is 3000) because the
  coarsest level is solved directly and sequentially (2.9 times faster on 128²); a smaller system
  AMGCL solves directly, in one iteration (coarsened to 500, the 1-D film of 2000 cells took 80
  iterations a solve and ran 6.5 times slower on one thread).
* **The concentration equation** with `ILU` tries a light preconditioner first (Jacobi with eigen,
  SPAI0 with amgcl), at most 300 iterations; a try that does not converge is dropped and the solve
  starts again from the same initial guess with ILUT or AMG. The same holds for the block ILUT. After
  a failed try the next 16 are skipped, doubling up to 1024, back to 16 after a try that converged.
  The iterations of a dropped try are counted in `JobInfo.txt`.
* **The tolerance of each Newton pass** of the film height follows the progress of the iteration
  (Eisenstat-Walker, choice 2): η_k = 0.9 (|F_k| / |F_k-1|)², safeguarded by 0.9 η_k-1² while that
  is above 0.1, between `Epsilon` and 1e-4; the first pass of a step at 1e-4. `Epsilon` is therefore
  the tightest tolerance a pass is solved to. On the cases measured the final film moved by 1e-15
  relative and the linear iterations fell by 24 to 38 %.

`YUZEN_OPT` turns each of these off (`-bilu`, `-amgcoarse`, `-concpc`, `-ew`; `none` all of them),
to compare with the algorithm it replaced.

An in-house BiCGStab with ILU(0) was measured during development and left out of the release: on
`leveling_1d_reg` it needed 3676 linear iterations where `eigen` needs 21, and its one property,
bit-identity under a change of thread count, is a property of `eigen` as well.

### 7.2 Which one to choose

Measured on one node with two Intel Xeon CPU Max 9480 (112 cores, 8 NUMA nodes), GCC 11 with LLVM's
`libomp` (the default build, §7.3), `Release` build with `-DYUZEN_WITH_AMGCL=ON`, threads bound to
cores 0..T−1, the shortest of three runs, up to the cores of one socket (§7.3).
**The absolute values depend on the machine, so what to read is the ratio between the columns.**

| Case | Threads | eigen | amgcl |
|---|---|---:|---:|
| 1-D `leveling_1d_reg` (Nx=300) | 1 | **0.58 s** | 0.64 s |
|  | 2 | **0.46 s** | 0.51 s |
|  | 4 | **0.38 s** | 0.44 s |
|  | 8 | **0.35 s** | 0.42 s |
|  | 16 | **0.36 s** | 0.42 s |
|  | 28 | **0.36 s** | 0.43 s |
|  | 56 | **0.37 s** | 0.44 s |
| 2-D `leveling_2d_reg` (60×60) | 1 | 19.8 s | **12.0 s** |
|  | 2 | 10.6 s | **7.5 s** |
|  | 4 | 6.5 s | **4.8 s** |
|  | 8 | 5.8 s | **3.6 s** |
|  | 16 | 5.5 s | **3.1 s** |
|  | 28 | 5.4 s | **2.9 s** |
|  | 56 | 5.3 s | **2.8 s** |

* **In 1-D eigen is 10 % faster on one thread and stays ahead; in 2-D amgcl is 1.6 times as fast on
  one thread and 1.9 times on 56.** These cases are small (300 and 3,600 cells). In 1-D both gain
  about 1.6 by eight threads and nothing after, and more threads cost nothing (hot teams, §7.3; 56
  is within 6 % of 8). In 2-D eigen gains 3.8 by 56 threads and amgcl 4.3, most of it by eight.
  §7.3 shows the larger films. Where the time of the 2-D case goes (`tools/profile_threads.py`):
  with eigen, building the block-Jacobi ILUT takes 55 % of the time-step loop on one thread and on
  eight, the BiCGStab iterations 13 %; the factorization is redone at every linear solve, and its
  blocks stop dividing the work at a few thousand rows. The profile of amgcl is spread out: on one
  thread AMG set-up 9 %, the iterations 27 %, the assembly of the film-height matrices 27 %.
* **If exact agreement between runs matters, use eigen.** AMGCL's internal reductions are not
  order-fixed, so changing the thread count can change the low-order bits. Measured
  (`leveling_1d_reg`, 1 thread against 4), the differences are 6.089e-11 in the physical quantities
  and 9.615e-2 in the error measures (a column that is rounding in that case); at 2 threads, and on
  `allphys_2d_reg` at 2 and 4, the bits happened not to move. Under the same conditions eigen gives
  exactly 0, by construction.

> **★ Changing the backend changes how the rounding falls ★**
> Results do not agree bit for bit; they agree only to within the convergence tolerance of the
> linear solver. When you make a regression comparison, fix the backend on both sides of it, and
> check the `Solver Info` section of each `JobInfo.txt`.

### 7.3 Shared-memory parallelism

**Yuzen has shared-memory parallelism (OpenMP) only and does not support distributed-memory
parallelism (MPI); no MPI version is planned.** One run is one process on one machine, using
several threads inside it, so the largest problem is bounded by the memory of one node, and a
cluster is useful as a way of running several cases at once, one process each.

The thread count comes from `OMP_NUM_THREADS`; there is no input key for it.

**The OpenMP runtime.** Every parallel loop costs a fork and a join, and a step makes thousands of
them. On a 56-core socket of a Xeon Max 9480 an empty parallel loop costs 8.0 µs on 56 threads with
GCC's `libgomp`, and 2.1 µs with LLVM's `libomp` or Intel's `libiomp5`
(`tools/scaling/omp_overhead.cpp`). So a GCC build links `libomp` instead of `libgomp` when it is
installed (`YUZEN_PREFER_LIBOMP`, default ON; `YUZEN_OPENMP_RUNTIME` names a runtime outright). With
`libomp` or `libiomp5`, a linear solve of n unknowns runs on min(T, ⌈n/500⌉) threads
(`util/Grain.h`; `YUZEN_GRAIN=lin=G` or `off`), and the solver sets `KMP_HOT_TEAMS_MODE=1` at
start-up unless it is set, which keeps resizing a team cheap. With `libgomp` resizing costs more than
it saves, and every region uses T. Measured on that socket, nine cases from 300 cells to 1024², both
backends: with `libomp` and grain 500 no case gets slower by more than 7 % up to 56 threads (with
`libgomp`, up to 3.6 times), and the fastest time is 1.15 times that of `libgomp` (geometric mean).
The run log and `JobInfo.txt` name the runtime (`OMP_RUNTIME`) and the grain (`YUZEN_GRAIN`).

**The thread limit.** Past one socket the reductions and the memory traffic cost more than the cores
return: every size up to 1024² ran slower on 112 threads than on 56. **A run uses at most the cores
of one socket** (from the Linux CPU topology; no limit where it cannot be read), and with `libgomp`
at most T = min(cores of one socket, max(8, ⌈cells/1000⌉)) (cells: Nx, Nx·Ny in 2-D, times Nz for
CZ), which came within 10 % of the fastest count there for every size measured; `util/ThreadCap.h`.
A request of 8 or fewer is never affected. **Asked for more with `OMP_NUM_THREADS`, the run stops**
before its result directory is made and says what to ask for; left to the OpenMP default, it takes
the limit. `YUZEN_THREAD_CAP` set to a number of cells per thread applies the size term with any
runtime, and `0` turns the limit off. `JobInfo.txt` records the count used (`OMP_MAX_THREADS`), the
count the runtime offered (`OMP_THREADS_REQUESTED`) and the cells per thread (`OMP_THREAD_CAP`).

* **Never ask for more threads than you have physical cores.** On a two-core machine,
  four threads on `leveling_1d_reg` cost a factor of three to six against one thread.
* **Cores that are busy count as missing.** The OpenMP runtime spins at the end of every parallel
  loop; with one of two cores busy, `concz_vs_wm` took 246 s on two threads against 5 s with
  `OMP_WAIT_POLICY=passive`. If you share the machine, use one thread per free core or set
  `OMP_WAIT_POLICY=passive` (the answer is the same either way); `ctest` sets it for the tests that
  run the solver.
* **On cores of your own, leave `OMP_WAIT_POLICY` unset.** `passive` puts a sleep and a wake-up into
  every parallel loop: an empty loop costs 25 µs on 8 threads and 258 µs on 112 with it, against
  1.7 and 24 µs unset (Xeon Max 9480, `tools/scaling/omp_overhead.cpp`).

What is parallelized is the assembly of the physical terms (the grid sweeps in `SetMatrixA` /
`SetJacobiMatrix`, about 60 % of a 1-D Newton iteration), the velocity calculation, the boundary
conditions, the reductions in `Vector` (in fixed blocks, §8.1), handing the matrix to the linear
solver (the gather into column order and the diagonal scaling, row by row), the block ILUT of eigen
and AMGCL's set-up and iterations, and the matrix-vector products. What remains sequential is the
direct solve of AMGCL's coarsest level (500 unknowns on a system of more than 3000) and the vector work of Eigen's BiCGStab.
Light loops run serially below 64 iterations a thread.

#### 7.3.1 Measured scaling

Wall-clock time of the time-step loop [s] on two Xeon Max 9480 sockets (112 physical cores, 8 NUMA
nodes, `OMP_PROC_BIND=close`), the thread limit off, the fastest of three repeats; measured at 1, 2,
4, 8, 14, 28, 56, 84 and 112 threads, and "best" is over all of them
(`tools/scaling/results/scaling-mc076.tsv`; the amgcl row of Nx = 2000 re-measured on another node
of the same kind after the fix to AMGCL's coarsest level, §7.2).

| Size | Backend | 1 | 8 | 28 | 56 | 112 | best |
|---|---|---:|---:|---:|---:|---:|---|
| 1-D, Nx=2000 | eigen | 9.745 | 4.079 | 3.651 | 4.329 | 5.466 | 2.7× at 28 |
| | amgcl | 10.59 | 4.056 | 3.66 | 4.32 | 5.426 | 2.9× at 28 |
| 1-D, Nx=8000 | eigen | 5.605 | 1.574 | 1.422 | 1.318 | 1.51 | 4.3× at 56 |
| | amgcl | 32.49 | 8.889 | 8.093 | 7.528 | 9.858 | 4.3× at 56 |
| 1-D, Nx=32000 | eigen | 8.363 | 2.909 | 2.18 | 2.022 | 2.2 | 4.1× at 56 |
| | amgcl | 19.98 | 3.978 | 2.657 | 2.678 | 3.787 | 7.5× at 28 |
| 1-D, Nx=131072 | eigen | 21.55 | 5.446 | 3.588 | 3.234 | 3.613 | 6.7× at 56 |
| | amgcl | 23.85 | 4.421 | 2.542 | 2.193 | 3.475 | 10.9× at 56 |
| 2-D, 128² | eigen | 27.52 | 5.386 | 3.506 | 3.275 | 3.328 | 8.4× at 56 |
| | amgcl | 27.07 | 5.354 | 3.585 | 3.379 | 3.613 | 8.0× at 56 |
| 2-D, 256² | eigen | 35.15 | 7.647 | 5.061 | 4.482 | 5.475 | 7.8× at 56 |
| | amgcl | 50.59 | 8.886 | 4.745 | 4.12 | 5.516 | 12.3× at 56 |
| 2-D, 512² | eigen | 58.1 | 12.0 | 7.144 | 6.29 | 9.164 | 9.2× at 56 |
| | amgcl | 95.18 | 14.23 | 5.542 | 4.226 | 9.138 | 22.5× at 56 |
| 2-D, 1024² | eigen | 147.9 | 26.92 | 15.13 | 13.09 | 19.01 | 11.3× at 56 |
| | amgcl | 170.7 | 24.94 | 11.06 | 12.03 | 9.209 | 18.8× at 84 |

* **Both backends scale, amgcl further**: at 56 threads eigen reaches 8–11× in 2-D and amgcl 8–22×,
  the larger factors on the larger films. Eigen's block-Jacobi ILUT is factorized in parallel; what
  limits it is the BiCGStab iterations (33–72 % of the loop), whose vector work is sequential.
* **Which is faster depends on the thread count.** In 1-D the two are level at Nx = 2000, where AMGCL
  solves the system directly, and eigen is ahead at every thread count from 8000 to 32000 (5.7× at 56 threads on Nx = 8000: AMGCL has little to coarsen and needs about 100
  linear iterations a solve where Eigen needs 12); at Nx = 131072 the two are within 11 % on one
  thread and amgcl is 1.5× faster at 56. In 2-D the two are level at 128²; from 256² up eigen is
  faster on one thread (1.2–1.6×) and amgcl from 8–28 threads up (1.1–1.5× at 56).
* **The useful thread count is 28–56, not 112.** Every case but amgcl on 1024² is at or past its
  optimum by 56 threads and slower at 112 than at 56; past one socket the reductions and the memory
  traffic cost more than the extra cores return.
* **Small problems do not parallelize.** At Nx = 2000 the best speed-up is 2.7× (eigen) and 2.9× (amgcl), both
  at 28 threads. Below a few thousand
  points, run one thread per case and run several cases side by side.

`tools/scaling/summarize.py` also prints an Amdahl fit $`f`$; where a case slows down past its
optimum it reports $`f = 0.000`$, which does not mean that nothing is parallel. Read the measured
speed-up and efficiency rather than $`f`$.

#### 7.3.2 Which backend to ask for

| Problem | Backend | Why |
|---|---|---|
| 1-D, up to ~10⁴–10⁵ cells | eigen | AMGCL has almost nothing to coarsen: 5.7× slower at Nx = 8000 |
| 1-D, ~10⁵ cells and more | either | within 11 % on one thread; amgcl 1.5× faster at 56 threads on Nx = 131072 |
| 2-D, a socket of threads | amgcl | 1.1–1.5× faster from 256² up, level at 128² |
| 2-D, eight threads or fewer | eigen | 1.2–1.6× faster from 256² up on one thread |

`auto` returns eigen for a 1-D film and amgcl for a 2-D one, so the default is the faster backend in
1-D below ~10⁵ cells and on a full socket in 2-D. In 2-D the two are within a factor of 1.6 up to one
socket; on a small machine and a 2-D film of 256² or more, set `Backend = "eigen"`. Where one overtakes
the other depends on the machine, so `auto` splits on the dimension alone; if it matters, measure both
(§7.3.3) and set `Backend` explicitly.

Changing the thread count does not change the answer with eigen (§8); with amgcl it moves the last
few digits. The thread count is a performance decision; the backend is the one with a correctness
dimension to it.

#### 7.3.3 Measuring it on your own machine (`tools/scaling/`)

`tools/scaling/` is a self-contained kit for this measurement. For details see
`tools/scaling/README.md`.

```sh
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release -DYUZEN_WITH_AMGCL=ON
cmake --build build -j
tools/scaling/scaling.sh --tier quick --repeats 1 --threads "1 4"   # smoke test
tools/scaling/scaling.sh --threads "1 2 4 8 16 32 64" --repeats 3   # full run
python3 tools/scaling/summarize.py tools/scaling/results/scaling-*.tsv
```

The measurement rests on the premise that **changing the thread count does not change the work**,
and the cases in `cases/scaling/` secure that in three ways: (1) the automatic time-step control is
switched off so the step is fixed, (2) the step count is fixed with `Numerics.MaxSteps` (`EndTime`
is set to 1e30), and (3) field output is switched off. (3) is not only to avoid timing the I/O:
**the output interval also caps the time step**. What is timed is `Time loop [s]` from
JobInfo.txt, that is, the wall-clock time of the time-step loop alone; `Elapsed Time` includes grid
construction and the initial output, so it is no guide on a large grid.

The results are appended one row per run to a TSV file, and combinations already measured are
skipped, so an interrupted sweep is restarted by repeating the same command. The head of the file
carries the CPU, the socket count, the core count, the NUMA node count, `OMP_PROC_BIND` and the
version of Yuzen, each prefixed with `#`.

**A small grid cannot fill a large machine, whatever you do.** The reductions in `Vector` are split
into a fixed number of blocks that does not depend on the thread count (§8.1), the number of blocks
being `min(256, N/512)`. Giving all p threads a block therefore requires `N >= 512p`, which at 112
threads is 57,344 points, or 131,072 points for a comfortable margin (the block count then reaches
its ceiling of 256). Conversely, **beyond 256 threads the reductions do not parallelize any
further**. This is the price paid by design for bit-level reproducibility.

### 7.4 The shared contract

* **The solution vector passed to `Solve()` is also the initial guess.** With Eigen it is handed to
  `solveWithGuess`. The return value is `0` for converged, `1` for the iteration cap reached, and
  `2` for divergence or breakdown.
* **A failed solve must have its solution vector thrown away.** If a NaN is left in the vector and
  passed to the next call, the solve fails with a NaN again no matter how far the time step is cut
  (`C_NewtonKantorovich::Solve` discards it; anything new that calls a backend directly must do the
  same).
* Diagonal scaling (`RowScaling`) divides each row by its own diagonal and so does not change the
  solution. The matrix itself is not modified; the scaling is applied to a mirror of the values and
  to a copy of the right-hand side.
* The preconditioner names map differently onto each implementation: under `amgcl`, `ILU`
  becomes smoothed-aggregation AMG with SPAI0, and `Jacobi` becomes damped Jacobi.
* The environment variable `YUZEN_ILUT_FILL` overrides the fill factor of Eigen's `IncompleteLUT`
  (default 10). Lowering it makes the preconditioner cheaper, but drop too much and the
  factorization breaks down and the solution comes back as NaN (a value of 1 breaks down even on a
  1-D problem).

## 8. Reproducibility

**This section is for developers who intend to modify the solver.**

Yuzen is written so that **changing the number of OpenMP threads does not change the result by a
single bit** (the `amgcl` backend excepted, §7.2). The tests assert this: `ctest` registers three
entries that run the same test binary in separate processes, with the environment's own thread
count, with `OMP_NUM_THREADS=1` and with `OMP_NUM_THREADS=4`. The property depends on the following
four things. Changing any of them for the sake of tidiness breaks bit-identity.

### 8.1 The deterministic reductions in `Vector`

`#pragma omp parallel for reduction(+:x)` decides how to split the partial sums, and in what order
to add them, from the number of threads. Floating-point addition is not associative, so a different
thread count gives different low-order bits in an inner product or a norm; those bits are then
amplified through the iteration count of BiCGStab and the acceptance test of the time-step control,
and surface as "one thread and two threads give different answers". The inner product
(`operator,`), `Sum()` and `Norm()` of `Vector` therefore split the array into **a fixed number of
blocks that does not depend on the thread count**, form a partial sum per block, and add the blocks
up in block order at the end. The number of blocks follows from the array length alone (at least
512 elements per block, at most 256 blocks), so the order of the additions is the same however many
threads are running.

> **★ Do not change the order of summation ★** Do not replace these with a bare
> `reduction(+:...)`. `tests/test_vector.cpp` pins this down.

`Max` / `Min` / `AbsMax` are exactly associative and exactly commutative, so a sequential loop is
fine for them.

### 8.2 The storage order of `CRS_Matrix`

Within a row the entries are stored **in insertion order**, and the first entry of every row is
always the diagonal. This is not a question of presentation: it fixes **the order of the
floating-point additions** in the matrix-vector product, in `ScalingDiag` and in `Residual`.

> **★ Do not sort the entries of a row into ascending order ★** Sorting changes how the sum rounds
> and bit-identity is lost.

Libraries that require ascending column indices within a row (Eigen, AMGCL) are given a **sorted
mirror**, kept separately from the matrix proper. The permutation is built once, at the moment the
structure is finalized; after that it is only an `O(nnz)` gather of the values.

Note that `CRS_Matrix` has two states: a construction phase (a variable-length array per row; the
rows are independent of one another, so insertion under OMP is safe) and a finalized phase (flat
CRS), and any read-side call finalizes it automatically. The BC flag is a `char` because
`std::vector<bool>` is a bit-packed proxy reference and is not thread-safe.

### 8.3 The registration order of the physical terms

The order in which `WellMixed_HeightEquationBase::RegisterPhysics()` lists the physical terms
decides the order of the `+=` operations into the matrix, that is, the column insertion order.

```
LaplacePressure → VerticalGravity → ParallelGravity
                → SolutoCapillarity → Centrifugal → DisjoiningPressure
```

> **★ Do not reorder the list; add a new term at the end ★**
> Reordering changes the storage order, and with it the rounding of the sums of §8.2. As long as a
> new term is appended at the end, no result of a calculation that does not use it changes by a
> single bit (a term that is not enabled does not even have its `SetMatrixA` called).

The velocity calculation (`CalVelocity`) is the one exception: each term writes into a velocity
field of its own there, so the order cannot affect the result.

### 8.4 The two direction orderings

The kernels carry two orderings of the directions and **use them in different places deliberately**
(`GridTraits.h`). `T::Dir(d)` is axis order (in 2D, Xm, Xp, Ym, Yp) and `T::StencilDir(k)` is
stencil order (in 2D, Ym, Xm, Xp, Yp). `SetMatrixA` sums each row in axis order, while
`SetJacobiMatrix` sums the diagonal and the cross terms in stencil order. The diagonal of the
capillary pressure term contains both in a single expression (`gBo*aB + gLo*aL + gRo*aR + gTo*aT`
is in stencil order, `sLo*bL + sRo*bR + sBo*bB + sTo*bT` is in axis order). In 1D the two orderings
coincide, so the distinction shows up only in 2D.

> **★ Follow the comment in the kernel about which ordering to use ★** Swapping one for the other
> breaks bit-identity. Do not "tidy up" by making both use the same ordering.

### 8.5 The remaining rules

* A preconditioner's triangular substitutions must run sequentially, because both the forward and
  the backward substitution carry a loop-carried dependency. Parallelizing them makes a thread read
  values another thread has not updated yet, which silently turns the preconditioner into a block
  Jacobi whose blocks are the thread decomposition, so the iteration count becomes thread-count
  dependent. The rule binds any preconditioner added to a backend of your own.
* `Vector` does not skip zero-initialization (`resizeFast()` does exactly what `resize()` does). Any
  path that can read an uninitialized value makes the result non-deterministic.
* Do not keep state in a function-local `static`. State that ought to belong to an instance becomes
  state belonging to the translation unit: several cases can no longer be run in one process, and
  counters bleed between separate instances (the thickness equation and the concentration
  equation).

## 9. The entries of the `Numerics` block

The defaults are those of `src/IO/Input.cpp`; the manual's chapter "Numerics and restart" is the
full reference (`ctest` entry `doc_keys` holds these tables to the parser).

### 9.1 `Numerics`

| Key | Default | Description |
|---|---|---|
| `MaxSteps` | `0` | Upper bound on the number of steps taken. At least 0; `0` sets no limit. Counted in 64 bits (up to 2^63−1); a restart counts on from the steps already taken |
| `InitialDt` | (required) | Initial time step [s], positive. Clamped down to the smallest output interval if it is larger |
| `EndTime` | `0` | End time [s], absolute (with a restart it must exceed `Restart.RestartTime`). Left at `0`, the run stops after a single step: always set it |
| `SteadyStateCriterion` | `1e-3` | Threshold of the steady-state test. Change rate per unit time $`\times`$ the span of the run. `0` or below disables the test |
| `SteadyStateSteps` | `10` | How many consecutive steps the criterion must hold (the test is `count > SteadyStateSteps`). Range `[1:1000000]` |

### 9.2 `Numerics.TimeStepControl`

| Key | Default | Description |
|---|---|---|
| `Enable` | `false` | Make the time step variable. When false, $`\Delta t`$ stays at `InitialDt` and a step that does not converge stops the run |
| `Criterion` | `0.01` | The local truncation error allowed per step (not the error itself, §5.1) |
| `ChangeFactor` | `0.5` | The factor multiplied in when shrinking (divided by when growing). Range `(0:1)`, both ends excluded |
| `KeepSteps` | `50` | How many steps the condition must keep holding before the step may grow. At least 1 |
| `MaxDtFactor` | `10000` | Upper bound on the step = this value $`\times`$ (minimum mesh width)$`^{4}`$, non-dimensional. Applied whether or not `Enable` is true; not applied by `CZ` |
| `MinDt` | `1e-15` | Lower bound on the step, non-dimensional. Requesting a smaller step aborts the run. Applied whether or not `Enable` is true; `CZ` uses a built-in `1e-15` |
| `ScaleFloor` | `0.0` | Floor on the weight of the error (as a fraction of the global maximum). Range `[0:1]` |

### 9.3 `Numerics.HeightEquation`

The `CZ` model solves the height equation semi-implicitly and ignores all five keys.

| Key | Default | Description |
|---|---|---|
| `NewtonMaxIterations` | `30` | Upper bound on the Newton outer iteration. Reaching it counts as a convergence failure, which shrinks the step when the step control is on and stops the run at a fixed step; so does a correction that stagnates far from convergence (§4.3) |
| `NewtonMinIterations` | `1` | Minimum number of iterations. Range `[1:100]` |
| `JacobianUpdateInterval` | `1` | How often the Jacobian is rebuilt. `1` rebuilds it every iteration. Range `[1:100]` |
| `NewtonEpsilon` | `1e-10` | Absolute convergence test (on the absolute value of the correction) |
| `NewtonEpsilonRel` | `1e-6` | Relative convergence test (as a ratio to the first correction). `0` or below disables it. Range `[0:1]` |

### 9.4 `Numerics.HeightEquation.LinearSolver` / `Numerics.ConcEquation.LinearSolver`

Both are required and have the same entries with the same defaults. `Numerics.ConcEquation` has no
keys of its own and is required even when no concentration equation is solved.

| Key | Default | Description |
|---|---|---|
| `Backend` | `"auto"` | Choice of implementation. `auto` / `eigen` / `amgcl`; `auto` is amgcl for a 2-D film on a build with AMGCL, otherwise eigen. Any other value is refused and stops the run (§7.1) |
| `PreconditionType` | `"ILU"` | Preconditioner. `None` / `Jacobi` / `ILU`; with `ILU` a cheaper one is tried first and this is the fall-back (§7.1) |
| `Epsilon` | `1e-13` | Threshold of the convergence test (relative residual). For the well-mixed film height, the tightest tolerance a Newton pass is solved to (§7.1) |
| `MaxIterations` | `5000` | Upper bound on the iteration. At least 1 |
| `RowScaling` | `true` | Apply diagonal scaling (does not change the solution) |

### 9.5 `Numerics.Convection`

Optional, at most once (§1.5). Both keys take one of `UTOPIA`, `KK`, `Koren`, `VanLeer`,
`Central4`, `QUICK`, `Upwind1`, `Upwind2` and `Central2`; any other name is refused with the list.
The schemes in use are written to the log and to `JobInfo.txt`.

| Key | Default | Description |
|---|---|---|
| `Concentration` | `"UTOPIA"` | The advection of $`c`$: in $`x`$ and $`y`$, and in `CZ` also in $`\zeta`$. One of `UTOPIA` / `KK` / `Koren` / `VanLeer` / `Central4` / `QUICK` / `Upwind1` / `Upwind2` / `Central2` |
| `Height` | `"Koren"` | The upwinded fluxes of the well-mixed film-height equation: the in-plane gravity, the centrifugal term and the Marangoni flux. Refused with `CZ`. One of `UTOPIA` / `KK` / `Koren` / `VanLeer` / `Central4` / `QUICK` / `Upwind1` / `Upwind2` / `Central2` |

The four diagnostic environment variables are `YUZEN_DEBUG_RESIDUAL` (§3), `YUZEN_DEBUG_COUPLING`
(§2.2), `YUZEN_DEBUG_FLOOR` (§4.1) and `YUZEN_ILUT_FILL` (§7.4).
