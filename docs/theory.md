# The theory behind Yuzen — governing equations and the meaning of the input parameters

This document explains which equations Yuzen solves and what each entry of the
input file refers to physically. It is a shorter companion to the user manual
(`docs/manual/`), whose chapters "The governing equations", "The physics
terms", "Wetting and the contact line" and "Concentration and boundary
conditions" give the full account, and whose reference chapters tabulate every
key with its type, default and range; where the two differ, the manual is
authoritative. Discretization, Newton iteration, linear solvers and time-step
control are in [`numerics.md`](numerics.md); comparisons against analytic
solutions and conservation laws are in [`validation.md`](validation.md).

---

## 1. Governing equations

### 1.1 The long-wave approximation

When the film height $`h`$ is sufficiently thin compared with the representative
in-plane length $`L`$, the Navier-Stokes equations collapse into one evolution
equation for $`h(x, y, t)`$. This is the long-wave (lubrication) approximation,
and it is the equation Yuzen solves. The four assumptions it makes are as follows.

* The aspect ratio is small, $`\varepsilon = h_0/L_0 \ll 1`$.
* Inertia is negligible. The velocity profile is parabolic across the thickness,
  with no slip at the substrate (or, with `[Physics.Slip]`, a Navier slip
  condition there; see below) and a shear stress at the free surface set solely
  by the surface-tension gradient.
* The pressure is uniform across the thickness. Capillary pressure, gravity and
  disjoining pressure all enter this single pressure alone.
* The interface slope is small, so the curvature may be linearized to $`\nabla^{2} h`$
  (the denominator of the exact mean curvature
  $`\nabla\cdot(\nabla h/\sqrt{1+\vert\nabla h\vert^{2}})`$ is set to 1).

**Navier slip at the substrate.** With `[Physics.Slip] Enable = true` the
no-slip condition $`u = 0`$ at $`z = s`$ is replaced by the Navier condition
$`u = \lambda\,\partial u/\partial z`$, with the slip length
$`\lambda`$ = `Physics.Slip.Length` [m]. The velocity profile is then the
parabola plus a uniform slip velocity, and integrating it across the film
gives

```math
\mathbf{Q} = -\frac{1}{\mu}\Big(\frac{h^{3}}{3} + \lambda h^{2}\Big)\nabla p
  \;+\; \frac{1}{\mu}\Big(\frac{h^{2}}{2} + \lambda h\Big)\nabla\sigma .
```

Every term driven by a pressure or a body force (capillary, disjoining,
vertical and in-plane gravity, centrifugal) has its mobility
$`h^{3}/3`$ replaced by $`h^{3}/3 + \lambda h^{2}`$, and the shear-driven
Marangoni term $`h^{2}/2`$ by $`h^{2}/2 + \lambda h`$; the viscosity
$`\mu(c)`$ divides both, as without slip. Everywhere below, read the mobility
that way when slip is on. Slip changes how fast the film moves, not where it
rests: an equilibrium (a uniform driving potential, §4.2) is the same with
and without it. Slip is implemented for the well-mixed model only; `CZ`
solves the velocity profile across the film itself and refuses the block.

### 1.2 The dimensional form

```math
\frac{\partial h}{\partial t} = -\nabla\cdot\mathbf{Q} \;-\; E(h, c) \;+\; S(x, t)
```

```math
\mathbf{Q} = \frac{\sigma h^{3}}{3\mu}\,\nabla\nabla^{2}(h+s)
  \;-\; \frac{\rho g_z h^{3}}{3\mu}\,\nabla(h+s)
  \;+\; \frac{h^{3}}{3\mu}\,\nabla\Pi(h)
  \;+\; \frac{h^{2}}{2\mu}\,\nabla\sigma
  \;+\; \frac{\rho g_\parallel h^{3}}{3\mu}
  \;+\; \frac{\rho\Omega^{2} h^{3}}{3\mu}\,r\,\mathbf{e}_r
```

$`h`$ is the film height and $`s(x, y)`$ is the height of the substrate, so the
free surface sits at $`h+s`$. Capillary pressure and vertical gravity are driven
by the shape of the free surface $`h+s`$, whereas disjoining pressure is a function
of the film height $`h`$ itself. Both $`\mu = \mu(c)`$ and $`\sigma = \sigma(c)`$ depend on the
concentration (§5.3). $`E`$ is the reduction of film height by evaporation and
$`S`$ is the supply from a nozzle. $`\mathbf{Q}`$ is the flux integrated across the
thickness. Nothing solves for a velocity: the in-plane velocity in the field
output is $`\mathbf{Q}/h`$, the depth-averaged velocity of the enabled terms,
and the `.vts` files of the well-mixed models also give it term by term, one
array per contribution except the disjoining pressure.

### 1.3 Concentration

The concentration $`c`$ of the solute (resin) is a mass fraction, so $`c = 1`$
represents the completely dry state. The default well-mixed model takes it to be
uniform across the thickness and evolves the solute $`hc`$ held in a column:

```math
\frac{\partial (hc)}{\partial t} + \nabla\cdot\left( c\,\mathbf{Q} \right)
  = \nabla\cdot\left( D(c)\,h\,\nabla c \right) \;+\; c_S\,S
```

The solute moves with the same flux $`\mathbf{Q}`$ as the liquid, diffuses along
the film, and arrives with the nozzle's liquid at the concentration $`c_S`$ that
`Physics.Source.ConcType` gives it: the concentration already there, $`c_S = c`$
(`"local"`, the default), or `Physics.Source.Conc` (`"supplied"`), mixed through
the thickness of the column (§3.8). Evaporation does not appear: the solvent leaves, the solute stays. This
is the equation the solver advances ([numerics.md](numerics.md) §1.7), so the
solute of a closed film is kept to rounding. Subtracting $`c`$ times the film
equation gives the same law for $`c`$,

```math
\frac{\partial c}{\partial t} + \mathbf{u}\cdot\nabla c
  = \frac{1}{h}\,\nabla\cdot\left( D(c)\,h\,\nabla c \right) \;+\; \frac{E c}{h}
    \;+\; \frac{S\,(c_S - c)}{h},
```

with $`\mathbf{u} = \mathbf{Q}/h`$ the mixed velocity; its $`E c/h`$ term is the effect
that "whatever the solvent leaves behind becomes more concentrated", stronger
where the film is thinner. The last term is 0 when the nozzle's liquid arrives at
the local concentration; liquid supplied at another $`c_S`$ pulls $`c`$ toward it
at the rate $`S/h`$, quickly on a thin film and slowly on a thick one. If the distribution across the thickness is
required, use the `CZ` model (§5.2).

### 1.4 Where the assumptions break down

* **Large slope.** The error of the linearized curvature, and of the
  $`1 - \cos\theta_e \approx \theta_e^{2}/2`$ used to calibrate the disjoining pressure, carries
  through directly. The error of the latter is 0.3 % at $`\theta_e`$ = 10°,
  1.0 % at 20°, 2.3 % at 30°, 5.3 % at 45° and 23 % at 90°.
  `Physics.DisjoiningPressure.ContactAngle` accepts 0 - 90°, but if the angle is
  to be used quantitatively it is better kept small.
* **When inertia matters.** At high rotation rates, or in a low-viscosity thin
  film, the velocity profile departs from the parabola.
* **When the in-plane length approaches the film height.** The neighborhood of a
  contact line is by construction a region of large slope, and the disjoining
  pressure is the device that allows it to be treated within the long-wave
  approximation (§4).

The assumption that the film is thin is concentrated in $`\mathrm{As} = (L_0/h_0)^{2}`$. A
configuration in which $`\mathrm{As}`$ falls to order 1 is outside the range of the
approximation.

## 2. Non-dimensionalization

### 2.1 Scales

Every input value is written dimensionally (m, s, Pa·s, ...), but the computation
is carried out in dimensionless variables.

| Symbol | Definition | Notes |
|---|---|---|
| $`L_0`$ | $`\sqrt{\sigma_0 / (\rho g)}`$ | Length (in-plane). `g = 9.80665 m/s²`, the capillary length |
| $`h_0`$ | `Initial.SurfaceShape.TotalHeight` | Height |
| $`T_0`$ | $`\mu_0 L_0^{4} / (\sigma_0 h_0^{3})`$ | Time |
| $`U_0`$ | $`L_0 / T_0`$ | Velocity |
| $`\varepsilon`$ | $`h_0 / L_0`$ | Aspect ratio |

$`\sigma_0`$ = `Fluid.SurfaceTension.BaseValue`, $`\rho`$ = `Fluid.Density.Value`,
$`\mu_0`$ = `Fluid.Viscosity.BaseValue`. Specifying `LengthScaleType = "custom"` makes $`L_0`$
the value of `LengthScale` instead of the capillary length.

The dimensionless film height is written $`H = h/h_0`$ and the in-plane coordinate
$`X = x/L_0`$. **The dimensionless slope $`H'`$ multiplied by $`\varepsilon`$ is the physical
slope**, and the whole discussion of contact angles (§4) rests on this conversion.
$`L_0`$, $`h_0`$, $`T_0`$ and $`U_0`$ are printed in the run log on standard
output under the heading `Scaling factors`; $`\varepsilon`$ is not printed, but
it is $`h_0/L_0`$ of the two that are.

### 2.2 Dimensionless groups

| Symbol | Definition | Physical meaning | Set by |
|---|---|---|---|
| $`\mathrm{As}`$ | $`(L_0/h_0)^{2}`$ | Inverse aspect ratio raised to the power 2; the strength of capillary pressure | `Initial.SurfaceShape.TotalHeight` and `L0` |
| $`\mathrm{Bo}`$ | $`\rho G_z L_0^{2} / \sigma_0`$ | Bond number; gravity normal to the substrate | `Physics.VerticalGravity.Acceleration` |
| $`\mathrm{Pgx}`$ | $`\rho G_x L_0^{2} / (\sigma_0\,\varepsilon)`$ | In-plane gravity (x); non-zero on a tilted substrate | `Physics.HorizontalGravity.Acceleration[0]` |
| $`\mathrm{Pgy}`$ | $`\rho G_y L_0^{2} / (\sigma_0\,\varepsilon)`$ | In-plane gravity (y) | `Physics.HorizontalGravity.Acceleration[1]` |
| $`\mathrm{Sp}`$ | $`\rho\Omega^{2} L_0^{4} / (\sigma_0 h_0)`$ | Spin number; centrifugal force | `Physics.Centrifugal.RotationRate` |
| $`\mathrm{Ma}`$ | $`(d\sigma/dc) / \sigma_0`$ | Marangoni number; how strongly surface tension depends on concentration | `Fluid.SurfaceTension.Coef` |
| $`\mathrm{Dp}`$ | $`\mathrm{As}\,(n-1)(m-1)\,\theta_e^{2} / (2\,H_*\,(n-m))`$ | Disjoining-pressure coefficient; fixed by the equilibrium contact angle | `Physics.DisjoiningPressure.*` and `PrecursorRatio` |
| $`H_*`$ | `PrecursorRatio` | Precursor thickness (dimensionless). $`\hat{\Pi}(H_*) = 0`$ | `Initial.SurfaceShape.PrecursorRatio` |
| $`\mathrm{Ev}`$ | $`E_0 T_0 / h_0`$ | Evaporation number; how fast evaporation thins the film | `Physics.Evaporation.Coef` |
| $`\mathrm{Pe}`$ | $`E_0 h_0 / D_0`$ | Peclet number; evaporation against diffusion | `Evaporation.Coef`, `Fluid.Diffusivity.BaseValue` |
| $`\mathrm{Df}`$ | $`D_0 T_0 / L_0^{2}`$ | Diffusion number; the speed of in-plane diffusion | `Fluid.Diffusivity.BaseValue` |
| $`\mathrm{Mr}`$ | $`M_\mathrm{resin} / M_\mathrm{solvent}`$ | Molar-mass ratio; used by Raoult's law | `Fluid.MolarMass.Solvent`, `.Resin` |
| $`\Lambda`$ | $`\lambda / h_0`$ | Dimensionless slip length; the mobility becomes $`H^{3}/3 + \Lambda H^{2}`$ and that of the Marangoni term $`H^{2}/2 + \Lambda H`$ (§1.1) | `Physics.Slip.Length` (0 unless `Physics.Slip.Enable`) |

$`G_z`$, $`G_x`$, $`G_y`$ are the components of gravitational acceleration given in the
input [m/s²], $`\Omega`$ is the rotation rate [rad/s], $`E_0`$ is the evaporation
coefficient [m/s] and $`D_0`$ is the reference diffusivity [m²/s].

$`\mathrm{Bo}`$, $`\mathrm{Pgx}`$, $`\mathrm{Pgy}`$, $`\mathrm{Sp}`$, $`\mathrm{Ma}`$, $`\mathrm{Dp}`$ **become 0** when the corresponding
term is disabled in the input, and the term is then not assembled at all. $`\mathrm{Ev}`$, $`\mathrm{Pe}`$, $`\mathrm{Df}`$ always carry a value,
and evaporation and diffusion are switched on and off by separate keys.
$`\mathrm{Ma} = 0`$ removes the Marangoni stress and, at the same time, removes the
concentration dependence of surface tension itself (§5.3).

If the capillary length is used as the length scale and `Acceleration = 9.80665` is set,
then $`\mathrm{Bo} = 1`$. `Acceleration` is the key to change when you want to move $`\mathrm{Bo}`$.

### 2.3 The dimensionless form

```math
\frac{\partial H}{\partial T} = -\nabla\cdot\mathbf{Q} \;-\; \mathrm{Ev}\,\hat{E}(H,C) \;+\; \hat{S}
```

```math
\mathbf{Q} = \frac{H^{3}}{3\hat{\mu}}\left[ \hat{\sigma}\,\nabla\nabla^{2}(H+S)
  \;-\; \mathrm{Bo}\,\nabla(H+S)
  \;+\; \mathrm{Dp}\,\nabla\hat{\Pi}(H)
  \;+\; \mathbf{Pg}
  \;+\; \mathrm{Sp}\,r\,\mathbf{e}_r \right]
  \;+\; \frac{\mathrm{As}\,\mathrm{Ma}}{2}\,\frac{H^{2}}{\hat{\mu}}\,\nabla C
```

```math
\hat{\Pi}(H) = (H_*/H)^{n} - (H_*/H)^{m}
```

and the solute balance

```math
\frac{\partial (HC)}{\partial T} + \nabla\cdot\left( C\,\mathbf{Q} \right)
  = \mathrm{Df}\,\nabla\cdot\left( \hat{D}(C)\,H\,\nabla C \right) + C_S\,\hat{S}
```

($`C_S = c_S`$, the concentration the supplied liquid arrives at; $`C`$ itself by default).

$`\hat{\mu}(C)`$, $`\hat{\sigma}(C)`$ and $`\hat{D}(C)`$ are the viscosity, surface tension and
diffusivity normalized to 1 at the reference concentration. The condition that gives $`\mathbf{Q} = 0`$, namely that the
contents of the square bracket be uniform, is the equilibrium condition, and it
is the starting point of §4. With Navier slip (§1.1) the prefactor $`H^{3}/3`$
becomes $`H^{3}/3 + \Lambda H^{2}`$ and $`H^{2}`$ of the Marangoni term
$`H^{2} + 2\Lambda H`$; the bracket, and with it the equilibrium condition,
is unchanged.

The bracket is the gradient of the driving potential
$`\varphi = \hat{\sigma}\nabla^{2}(H+S) + \mathrm{Dp}\,\hat{\Pi}(H) - \mathrm{Bo}\,(H+S)`$
(for $`\hat\sigma`$ uniform), and $`\varphi`$ is minus the variation of the free energy

```math
\mathcal{F} = \int\Big[\tfrac12\hat{\sigma}\,\vert\nabla(H+S)\vert^{2}
  + \tfrac{\mathrm{Bo}}{2}(H+S)^{2} + \mathrm{Dp}\,w(H)\Big]dA,\qquad
w(H) = -\int_{H_*}^{H}\hat{\Pi}(s)\,ds .
```

With the capillary pressure, the vertical gravity and the disjoining pressure
alone, a uniform concentration and closed or periodic edges,
$`d\mathcal{F}/dt = -\int \frac{M}{\hat\mu}\vert\nabla\varphi\vert^{2}dA \le 0`$
($`M = H^{3}/3`$, or $`H^{3}/3 + \Lambda H^{2}`$ with slip): the film runs downhill
in $`\mathcal{F}`$ to an equilibrium. `[Output.Energy]` writes $`\mathcal{F}`$, the
dissipation and their balance (`Energy.dat`). Film height and concentration live at cell centers
and velocities on cell faces.

## 3. The individual terms

Each term has its own independent sub-block under the `Physics` block and can be
switched on and off 1 at a time with `Enable`. **A sub-block left out of the
file is a term switched off**: it reads as an empty block with its switch off,
and every other key takes its default. For the four sub-blocks whose `Enable`
defaults to `true` — `LaplacePressure`, `SolutoCapillarity`, `Evaporation` and
`Source` — leaving the block out is therefore not the same as writing it
empty: an empty `[Physics.Evaporation]` switches evaporation on, a missing one
leaves it off. A missing `ConcField` means `Model = "None"`.

### 3.1 Capillary (Laplace) pressure

```math
\mathbf{Q}_\mathrm{cap} = \frac{\hat{\sigma}(C)\,H^{3}}{3\hat{\mu}(C)}\,\nabla\nabla^{2}(H+S)
```

The effect of surface tension acting to flatten the free surface. The leveling
of a film, a droplet relaxing toward a spherical cap, and liquid collecting on a
textured substrate are all this term. It contains a 4th-order spatial derivative,
and it is because of this term that the film-height equation is 4th order. The
input key is `Physics.LaplacePressure.Enable`. The curvature is measured on the
free surface $`H+S`$, not on the film height $`H`$. The mobility contains the
concentration dependence of the surface tension as well as that of the viscosity,
so if drying changes $`\sigma`$ the rate of leveling changes with it.

### 3.2 Disjoining pressure

```math
\mathbf{Q}_\mathrm{dis} = \mathrm{Dp}\,\frac{H^{3}}{3\hat{\mu}}\,\nabla\hat{\Pi}(H), \quad
\hat{\Pi} = (H_*/H)^{n} - (H_*/H)^{m} \quad (n > m > 1)
```

Once a film becomes thin, the interaction between the substrate and the free
surface starts to take effect. Expressing this as a pressure gives the disjoining
pressure, and it is this term that carries partial wetting (a finite contact
angle) and the motion of the contact line. The 1st term is repulsive (it tries
to keep the film thick) and the 2nd is attractive. $`\hat{\Pi} = 0`$ at $`H = H_*`$, and
this film height is what remains outside a droplet as the precursor film. The
input keys are `Enable` / `ContactAngle` / `ExponentN` / `ExponentM` in
`Physics.DisjoiningPressure`, together with
`Initial.SurfaceShape.PrecursorRatio`, which sets the precursor thickness.
The exponents default to $`(n, m) = (3, 2)`$, whether the block is absent or
written without them. Enabling this term with `PrecursorRatio = 0` leaves
$`H_*`$ undefined and terminates the run before any step, with
`Physics.DisjoiningPressure requires Initial.SurfaceShape.PrecursorRatio > 0.`
in `FatalErrors.txt`. It carries the strong non-linearity $`\Pi \sim H^{-n}`$ and is the stiffest
term in the model. The details are in §4.

### 3.3 Vertical gravity and the substrate shape

```math
\mathbf{Q}_\mathrm{vg} = -\frac{\mathrm{Bo}}{3}\,\frac{H^{3}}{\hat{\mu}}\,\nabla(H+S)
```

Gravity normal to the substrate. It drives liquid from the high parts of the free
surface to the low ones. On a horizontal substrate it acts to level bumps and
contributes to flattening alongside capillary pressure. The input keys are
`Physics.VerticalGravity.Enable` and `Acceleration` [m/s²]. If the substrate is used
tilted by an angle $`\alpha`$, give `Acceleration` = $`g\cos\alpha`$.

The substrate shape $`S(x, y)`$ is given in the `Substrate` block. List as many
`RectBump`, `EllipseBump`, `PeriodicWave`, `FlatPlane`, `MexicanHat` and
`SmoothConvex` blocks as required; they are superposed in the order of `Order`
into 1 single shape. The substrate does not change with time. $`S`$ enters the
capillary pressure and the vertical gravity in the combination $`H+S`$, so unless
one of these two terms is enabled the substrate shape does nothing.

### 3.4 In-plane gravity

```math
\mathbf{Q}_\mathrm{hg} = \frac{\mathbf{Pg}}{3}\,\frac{H^{3}}{\hat{\mu}}
```

Flow down an inclined substrate. What distinguishes it from vertical gravity is
that the driving field is a constant vector independent of position; since a
thicker film flows faster, the thick parts run ahead and form a front. The input
keys are `Physics.HorizontalGravity.Enable` and `Acceleration = (Gx, Gy)` [m/s²]. Give
the in-plane component $`g\sin\alpha`$ corresponding to the tilt angle $`\alpha`$ (`4.903`
at 30°). The sign sets the downhill direction.

### 3.5 Centrifugal force

```math
\mathbf{Q}_\mathrm{cf} = \frac{\mathrm{Sp}}{3}\,\frac{H^{3}}{\hat{\mu}}\,r\,\mathbf{e}_r
```

The effect of the liquid film being thrown outwards when the substrate rotates;
the main driving force of spin coating. The driving field is the distance $`r`$ from
the rotation axis, so it is stronger toward the rim. The input keys are
`Physics.Centrifugal.Enable`, `RotationRate` [rad/s] and `RotationAxis` (the
coordinates of the rotation center [m]). In a 1-dimensional (`xh`) calculation,
the material that would be thrown into the unsolved y direction and leave is also
treated as outflow, $`\tfrac{\mathrm{Sp}}{3}H^3/\hat{\mu}`$, in both models (`CZ`
also carries it into its vertical velocity); `History.dat` counts it, and the
solute it carries, in `OutflowedVol` and `OutflowedResin`, so that `VolErr`
balances. A rotating film is not periodic (the arm $`r`$ jumps by the length of
the domain across a seam), so a `Periodic` edge is refused with this term,
except in a run with `Numerics.ManufacturedSolution`.

### 3.6 Solutal capillarity (Marangoni)

```math
\mathbf{Q}_\mathrm{sc} = \frac{\mathrm{As}\,\mathrm{Ma}}{2}\,\frac{H^{2}}{\hat{\mu}}\,\nabla C
```

A surface stress arising because surface tension depends on concentration. Where
drying is faster the concentration rises, the surface tension changes, and the
resulting gradient drags the liquid. It is proportional to $`H^{2}`$ rather than
$`H^{3}`$ because it is a shear applied at the surface rather than a body force. It
is the origin of drying mottle and flow streaks. The input keys are
`Physics.SolutoCapillarity.Enable` and `Fluid.SurfaceTension.Coef` ($`d\sigma/dc`$ [N/m]).
The sign of `Coef` decides which way the liquid is pulled. `Enable = false` forces
$`\mathrm{Ma} = 0`$, so as this term disappears the surface tension also ceases to depend on
concentration.

### 3.7 Evaporation

| Equation | Evaporation term |
|---|---|
| Height equation | $`-\mathrm{Ev}\,\hat{E}(H, C)`$ |
| Solute $`HC`$ (the equation solved) | none |
| Concentration $`C`$ (the advective form) | $`+\mathrm{Ev}\,\hat{E}(H, C)\,C / H`$ |

Solvent leaves, the film thins, and the solute left behind becomes more
concentrated. Written for $`C`$, the term is divided by $`H`$, so
concentration rises fastest where the film is thinnest. The evaporation rate is
the product of a "part that depends on concentration" and a "part that depends on
position".

```math
\hat{E} = \text{Coef}\cdot f_\mathrm{conc}(C)\cdot f_\mathrm{space}(X)\cdot f_h(H)
```

* `ConcFunc` — `"Raoult"` is Raoult's law, proportional to the mole fraction of solvent (computed from the
  ratio $`\mathrm{Mr}`$ of `Fluid.MolarMass`). `"PowerN"` is
  $`\max(1 - C - \mathrm{ConcR}, 0)^{\mathrm{PowerN}}`$, which never goes
  negative.
* `SpatialFunc` — `"uniform"` is uniform in the plane. `"negPow"` is a radial
  distribution appropriate to a rotating substrate,
  $`\min(k\,(1-r/R')^{-n}, 10)`$ with $`k = \tfrac12(1-n)(2-n)`$ and
  $`R' = 1.001\,R`$; for $`0 < n < 1`$ it diverges toward the rim of the pool as a
  negative power of the distance from it, $`n = 0`$ is uniform, $`n < 0`$ vanishes at the rim. Over 5 cells it
  is blended into a quarter cosine that takes over the curve's value and slope,
  peaks at `NegPowRmax` (in meters, measured from the origin of the grid
  coordinates) and is constant beyond it. `NegPowExponent` sets $`n`$ (less than 1).
* $`f_h`$ — a correction that suppresses evaporation from the precursor film. It is
  applied only when `PrecursorRatio > 0`, and ramps up linearly to be 0 at
  $`H = H_*`$ and 1 at $`H = 20 H_*`$; it is 0 below $`H_*`$ as well, so a
  film thinner than the precursor neither evaporates nor takes up solvent. Its purpose is to stop the precursor film
  drying out and disappearing, so it is always active in a calculation that
  contains a contact line.

**No step removes more solvent than the film holds**: $`E\Delta t \le h(1-c)`$, at
either end of the step and in both equations, beside the cap at
$`h/(2\Delta t)`$. A film therefore dries toward $`c = 1`$ without passing it
— with a rate that does not depend on $`c`$, a third of the remaining solvent is
left after each step — and no solute is ever removed.

`Coef` has the dimensions of a velocity [m/s] and is the source of both $`\mathrm{Ev}`$ and
$`\mathrm{Pe}`$. The condition under which the run is judged dry and stopped is set by
`EndConcType` (`"mean"`, `"max"` or `"min"`; any other string is refused) and
`EndConcValue` (default 0.95). The dry-out test runs whether or not evaporation
is enabled, and also when the block is left out.

### 3.8 The source (nozzle)

The supply term $`\hat{S}`$ added to the height equation. It is modeled as liquid
dispensed from a nozzle mounted on the end of a rotating arm. The input keys are
in the `Physics.Source` block.

* `Radius` — the radius of the nozzle mouth [m]. A cell receives supply in
  proportion to the fraction of it covered by that circle (by that interval in 1
  dimension). `SubRegions` is the number of sub-divisions per cell used to
  measure that fraction; a sub-division counts when its center is inside.
* `InletRate` — the steady-state supply rate, given with the dimensions of a rate
  of thickness increase [m/s]. `TimeConstant` is the rise-time constant $`\tau`$ [s],
  and the rate grows as $`1 - \exp(-t/\tau)`$.
* `ArmCenter`, `ArmLength`, `RotationRate`, `StartAngle` — the center of the arm
  [m], its length [m], the rotation rate [rad/s] and the initial angle [rad].
  `ArmLength = 0` gives a nozzle that does not move. `RampTime` / `PlateauTime` are
  the time applied to the ramp-up and ramp-down of the rotation rate, and the time
  it holds constant [s]: the rate is $`\omega f(t)`$, with $`f`$ rising linearly
  from 0 to 1 over $`T_s`$, holding 1 for $`T_m`$ and falling back to 0 over
  another $`T_s`$. The arm's angle is the integral of that rate,
  $`\theta_0 + \omega\int_0^t f\,dt'`$, so it stops at
  $`\theta_0 + \omega(T_s + T_m)`$ and stays there.

`ConcType` says what the supplied liquid brings to the solute balance, as
$`c_S S`$ (well-mixed model; `CZ` ignores the nozzle):

* `"local"` (the default) — the liquid has the concentration already present
  where it lands, $`c_S = c`$: it thickens the film without changing $`c`$, and
  `Conc` is not used.
* `"supplied"` — the liquid has the concentration `Conc` (a mass fraction in
  [0, 1]; a value outside is refused), mixed through the thickness of the column
  it lands on. On a film that does not flow,
  $`c(t) = (h_0 c_0 + c_S \int_0^t S\,dt')/(h_0 + \int_0^t S\,dt')`$, so a thin
  film under the nozzle takes on $`c_S`$ quickly and a thick one slowly.

`History.dat` counts the solute the liquid brings in `SuppliedResin`, at the
concentration the mode gives it.

## 4. Wetting and the contact line

### 4.1 The disjoining-pressure model and how Dp is fixed

```math
\Pi(h) = B \left[ (h_*/h)^{n} - (h_*/h)^{m} \right] \quad (n > m > 1)
```

```math
B = \frac{(n-1)(m-1)}{h_*\,(n-m)}\,\sigma\,(1 - \cos\theta_e)
```

$`h_*`$ is the precursor thickness and $`\theta_e`$ is the target equilibrium contact angle.
The default exponents are $`(n, m) = (3, 2)`$. The coefficient $`B`$ is fixed by the
equilibrium condition "difference in surface energy = integral of the disjoining
pressure", namely $`\sigma\,(1 - \cos\theta_e) = B h_* (n-m) / ((n-1)(m-1))`$. In the long-wave
approximation $`\theta_e`$ is small, so $`1 - \cos\theta_e \approx \theta_e^{2}/2`$ is substituted, and dividing
by the pressure scale $`\sigma h_0 / L_0^{2}`$ gives the coefficient the code actually
builds,

```math
\mathrm{Dp} = \frac{\mathrm{As}\,(n-1)(m-1)}{2\,H_*\,(n-m)}\;\theta_e^{2}
```

Changing `ContactAngle` moves $`\mathrm{Dp}`$ as $`\theta_e^{2}`$. Changing `ExponentN` /
`ExponentM` or `PrecursorRatio` is compensated on the side of $`\mathrm{Dp}`$, so **the
target contact angle does not change**. What changes is the thickness of the
precursor film and the width of the transition region.

That this expression is consistent with the contact angle can be verified
analytically. In the limit in which the precursor extends to infinity
($`K = 0`$, $`H_\mathrm{f} = H_*`$),

```math
\int_{H_*}^{\infty} \hat{\Pi}\,dH = \frac{H_*(m-n)}{(n-1)(m-1)}
```

```math
(H')^{2} \to \frac{2\,\mathrm{Dp}\,H_*\,(n-m)}{(n-1)(m-1)} = \mathrm{As}\,\theta_e^{2}
```

so the physical slope, the dimensionless slope $`H'`$ multiplied by $`\varepsilon`$, is exactly
$`\theta_e`$.

### 4.2 The equilibrium droplet

The equilibrium condition is zero flux, that is, uniform pressure; multiplying by
$`H'`$ and integrating once gives a first integral.

```math
H'' + \mathrm{Dp}\,\hat{\Pi}(H) = K
```

($`K`$ is a constant = the pressure)

```math
\frac{1}{2} H'^{2} = K (H - H_\mathrm{f}) - \mathrm{Dp} \left[ W(H) - W(H_\mathrm{f}) \right], \quad
W(H) = \int_{H_*}^{H} \hat{\Pi}\,ds
```

$`H_\mathrm{f}`$ is the flat film surrounding the droplet, where $`H' = H'' = 0`$ and hence
$`K = \mathrm{Dp}\,\hat{\Pi}(H_\mathrm{f})`$. At the summit $`H_\mathrm{max}`$ too $`H' = 0`$, so
$`K (H_\mathrm{max} - H_\mathrm{f}) - \mathrm{Dp} \left[ W(H_\mathrm{max}) - W(H_\mathrm{f}) \right] = 0`$ ties $`H_\mathrm{f}`$ to $`H_\mathrm{max}`$. In other words
**equilibrium droplets are a 1-parameter family**: fix $`H_\mathrm{max}`$ and the shape can
be built from 1 single ordinary differential equation.
`tools/make_equilibrium_1d.py` constructs this shape and writes it out in the
ASCII restart format.

### 4.3 The apparent contact angle and the precursor bias

The equilibrium droplet of this model, which has a finite precursor, does not show
$`\theta_e`$ unless $`H_\mathrm{max}/H_*`$ is sufficiently large. Outside the droplet (the region
where disjoining pressure has no effect) the profile is the parabola $`H'' = K`$,
and the slope obtained by extrapolating it down to $`H = 0`$ corresponds to the
contact angle measured experimentally. It can be written in closed form.

```math
\theta_\mathrm{cap} = \varepsilon\,\sqrt{-2 K H_\mathrm{max}} = \theta_e \left( 1 - \frac{H_*}{2 H_\mathrm{max}} \right)
```

The values for $`\theta_e = 10^\circ`$ and $`(n, m) = (3, 2)`$.

| $`H_\mathrm{max}/H_*`$ | apparent angle [deg] | max. slope of the profile [deg] |
|---:|---:|---:|
| 2 | 7.50 | 0.63 |
| 5 | 9.00 | 3.60 |
| 10 | 9.50 | 5.42 |
| 20 | 9.75 | 6.76 |
| 50 | 9.90 | 7.96 |
| 100 | 9.95 | 8.57 |
| 1000 | 9.995 | 9.55 |

At $`H_\mathrm{max}/H_* = 2`$ the apparent angle is 7.5°, and the maximum slope of the profile
is no more than 0.63°. **The 10° that was specified is simply not present in the
model under these conditions.** A configuration of "a 10 µm droplet on top of a
2 µm precursor film" does not stand up as a calculation for measuring a contact
angle. The convergence of the "maximum slope of the profile" is much slower than
that of the apparent angle, so measuring the angle from the slope requires a
thicker droplet still. This table can be reproduced with
`tools/make_equilibrium_1d.py --out /dev/null --sweep` (the sweep writes no
file, but `--out` is a required argument).

### 4.4 Rules of thumb for setting up a case

The length of the transition region over which the profile changes from the
precursor film to the body of the droplet is

```math
\ell = h_* / \theta_e
```

($`\theta_e`$ in radians)

For $`\theta_e = 10^\circ`$ and $`h_* = 0.2\ \mu\mathrm{m}`$ this gives $`\ell = 1.15\ \mu\mathrm{m}`$. 2 requirements follow
from it.

* **To resolve the contact line, aim for $`\Delta x \le \ell/5`$.** If this is not resolved,
  neither the speed at which the contact line moves nor the angle comes out right.
* **To obtain the angle you asked for, aim for $`H_\mathrm{max}/h_* \ge 50`$.** As the table in
  §4.3 shows, this ratio keeps the error in the apparent angle within 1 %. A ratio
  of 10 leaves you 5 % out, and a ratio of 5 leaves you 10 % out.

These 2 together fix the number of grid points. For instance $`h_0 = 10\ \mu\mathrm{m}`$, $`h_* = 0.2\ \mu\mathrm{m}`$
(`PrecursorRatio = 0.02`), $`\theta_e = 10^\circ`$, a domain of 600 µm and `Nx = 2000` ($`\Delta x = 0.3\ \mu\mathrm{m}`$)
is a configuration that sits right at the **borderline**. $`H_\mathrm{max}/h_* = 50`$ meets the
2nd rule exactly, but $`\Delta x = 0.3\ \mu\mathrm{m}`$ is slightly coarser than $`\ell/5 = 0.23\ \mu\mathrm{m}`$, so the
1st is not met (it amounts to $`\ell/3.8`$). Whether this is sufficient is answered
by the measurements in [`validation.md`](validation.md) §6.2: on this mesh the
end of spreading differs from an independently constructed equilibrium shape by
0.02 %, and so it does at the four times coarser $`\Delta x = 1.2\ \mu\mathrm{m}`$
(comparable to $`\ell`$), where the spreading itself follows the fine mesh's to
0.06 %. That case is lenient (a 10° droplet whose foot is far from any other
feature); $`\Delta x \le \ell/5`$ is no more than a rule of thumb, and whether a
coarser mesh suffices is decided by refining. Lowering
$`h_*`$ to gain $`H_\mathrm{max}/h_*`$ also shortens $`\ell`$ and therefore the required $`\Delta x`$; this
trade-off is what sets the cost of a calculation containing a contact line. If
evaporation is used as well, note that evaporation from the precursor film is
suppressed by the correction of §3.7, that is, $`h_*`$ enters the drying behavior
too.

## 5. The treatment of concentration

### 5.1 well-mixed (`WM`)

`Physics.ConcField.Model = "WM"` assumes the concentration is uniform across the
thickness. There is 1 unknown $`c(x, y, t)`$ per point and no mesh in the
thickness direction at all. It is sound when diffusion across the thickness is
much faster than evaporation, that is, when $`\mathrm{Pe} = E_0 h_0 / D_0`$ is small. When $`\mathrm{Pe}`$ is
large, only the region near the free surface becomes concentrated, and in reality
a hard layer (a skin) forms at the surface which slows all subsequent drying,
which well-mixed cannot represent. Because it under-predicts the surface
concentration, both the rise in viscosity and the Marangoni stress come out too
weak. The equation is solved as the balance of the solute $`hc`$ of each cell,
the solute crossing a face being the concentration there times the very volume
flux the film-height equation puts through it: the solute of a closed film is kept
to rounding, with evaporation or without, a uniform $`c`$ stays uniform, and $`c`$
is kept within [0, 1] by moving solute between cells, never by cutting it
([numerics.md](numerics.md) §1.5, §1.7). `EnableConvection` / `EnableDiffusion`
switch off the convective and in-plane diffusive terms individually.

### 5.2 `CZ` — resolving the thickness direction

`Physics.ConcField.Model = "CZ"` solves for the concentration with the thickness
direction resolved: $`c(x, z, t)`$ with `Grid.Dimension = "xh"`, and
$`c(x, y, z, t)`$ with `"xyh"` (the pseudo-3D case; `cases/concz_2d.toml` is a
worked example). The film is mapped onto the fixed coordinate
$`\zeta \in [0, 1]`$, so as the film thins the $`\zeta`$ mesh contracts with it
and the number of cells through the thickness, `Grid.Nz`, stays constant.

* The concentration equation is solved **in conservative form for $`hc`$** on the
  $`\zeta`$ grid: each cell exchanges solute through face fluxes, and the flux
  through its $`\zeta`$ faces follows from the discrete conservation of the grid
  with the film's step, so a uniform concentration stays uniform and the solute
  of a closed film is kept to rounding.
* **No solute crosses the substrate or the free surface**: the total flux there —
  the advection by the evaporating solvent and the whole diffusive flux,
  including its part along the slope of the surface or the substrate — is zero.
  The solvent that leaves the top cell without its solute is what builds up the
  concentrated layer near the surface. The evaporation rate uses **the
  free-surface concentration**, not the depth average, which is what lets the
  surface layer slow its own drying.
* Because the viscosity varies across the thickness, the in-plane velocity cannot
  be written down analytically, and the velocity profile is built by integrating
  $`1/\mu`$ and $`z/\mu`$ numerically across the thickness. The film height is
  advanced in the same conservative form as the well-mixed model, with the depth
  integral of that profile as its flux (the kinematic condition at the free
  surface): the mobility carries the mean of $`1/\mu(c(\zeta))`$ weighted by
  $`3(1-\zeta)^2`$, the Marangoni flux the mean weighted by $`2(1-\zeta)`$, and
  the surface tension is that of the free-surface concentration.

The costs are that the number of concentration unknowns grows from $`N_x`$ to
$`N_x \times N_z`$ and that the physics terms available are limited. **`CZ` has no
disjoining pressure and no source term**: `Physics.DisjoiningPressure` and the
whole of `Physics.Source` are ignored without a diagnostic. What is available is
capillary pressure, solutal capillarity, vertical gravity, in-plane gravity,
centrifugal force and evaporation, and every calculation involving a contact
line is done with `WM`. `CZ` also implements only `Symmetry` and `Periodic`
edges (§6).

### 5.3 Concentration dependence of the material properties

The `Fluid` block gives the value at the reference concentration `BaseConc` and
the law by which it varies away from it.

| Expression | Notes |
|---|---|
| $`\sigma(c) = \text{BaseValue} + \text{Coef}\cdot(c - \text{BaseConc})`$ | Surface tension, linear |
| $`\mu(c) = \text{BaseValue}\cdot\exp\left( \text{Exponent}\cdot(c - \text{BaseConc}) \right)`$ | Viscosity, exponential |
| $`D(c) = \text{BaseValue}\cdot\exp\left( \text{Exponent}\cdot(c - \text{BaseConc}) \right)`$ | Diffusivity, exponential |

That the viscosity rises exponentially as drying proceeds and the concentration
rises is the strongest non-linearity in this solver. The `Exponent` of the
viscosity can be around 10, and at that value a change in concentration of just
0.1 multiplies the viscosity by e. This is the effect by which the film solidifies
and the flow stops. The slope `Coef` of the surface tension is the source of $`\mathrm{Ma}`$
and at the same time enters the mobility of the capillary pressure as $`\sigma(c)`$. The
viscosity enters every flow term as $`\mu(c)`$, so the concentration field is coupled
to the film-height equation in both directions.

## 6. What the boundary conditions mean

A `Type` is specified for each of `XminBC` / `XmaxBC` / `YminBC` / `YmaxBC` in
the `Boundary` block. The height equation is 4th order, so 1 boundary needs
two conditions, and each type names a pair; the concentration equation is 2nd
order and takes one condition from the same name. What follows is the physical
meaning; the conditions each type imposes on the two ghost layers are given in
the manual's chapter "Boundary conditions and output", and how the ghost rows
enter the matrix in [`numerics.md`](numerics.md) §1.3.

* **`Symmetry`** — a plane of symmetry. Nothing crosses the boundary, and both the
  shape and the concentration are mirrored. It serves both for putting the center
  of a droplet or a film on the boundary to halve the cost, and for use as a closed
  wall. The default.
* **`Periodic`** — periodic. The edge is regarded as joined to the opposite one.
  Use it for a film on a periodically structured substrate, or to compute just 1
  wavelength of a pattern that continues indefinitely. It must be specified at
  both ends of a pair (specifying it on one side only, in 1-D or 2-D, is refused
  before the run starts), and it is refused with the centrifugal term (§3.5).
* **`Pinning`** — pinning. The film thickness at the boundary is held at zero and
  flow across the boundary is stopped, so the initial film must reach zero there
  too: before the first step the solver extrapolates the film to the wall,
  reports the value in the run log, and stops the run if more than a tenth of
  the thickest film is left there (`Initial.SurfaceShape.EdgeBump` is the shape
  made for such an edge). On a restart the film at the wall is taken as it is,
  and the same line is only a warning. It represents a situation in which
  the edge of the liquid is arrested by a step or by contamination, without
  resolving a contact line. It is the condition for computing "a film whose edge
  does not move" without using the disjoining pressure. The concentration takes a
  zero-gradient condition.
* **`FreshInlet`** — inflow. There is no key for the inlet state: the
  concentration at the boundary is held at the mean of the two points either
  side of the face in the field the boundary condition was built from, that is,
  at its initial value for the whole run; the film height takes the same mean
  from the current field, recomputed at every assembly, so the composition of
  what arrives is fixed and its thickness follows the solution. It represents
  liquid of constant composition arriving at the upstream end of a coating
  operation.
* **`Outlet`** — free outflow. The high-order derivatives are set to zero so that
  the liquid leaves without changing shape. Use it where the domain was merely cut
  for convenience and you have no information about what lies beyond. The volume
  that leaves is reported in `History.dat` (`OutflowedVol`).

Two combinations are refused or not implemented. The `CZ` model implements
`Symmetry` and `Periodic` edges only, and the input check refuses `Pinning`,
`FreshInlet` and `Outlet` with it, naming the edge. In 2-D, a corner where
`Pinning` meets `FreshInlet` or `Outlet` is not implemented and aborts the run
(`not yet implemented` in `FatalErrors.txt`).

The `Substrate` shape is mirrored at every edge except a `Periodic` one, where
it is continued from the other end, whatever condition is specified for the
film height and the concentration; place a wall where the substrate is flat if
you can. The choice of boundary has a
strong effect on the result. In a calculation where a droplet spreads and reaches
the boundary, `Symmetry` is the same as having a wall there, and gives a different
answer from `Outlet`.
