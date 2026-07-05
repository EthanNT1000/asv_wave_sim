# Hydrodynamics Physics Algorithm

**Source Files**: [`gz-waves/src/Physics.cc`](https://github.com/EthanNT1000/asv_wave_sim/blob/ament_environment_hooks/gz-waves/src/Physics.cc), [`gz-waves/src/systems/hydrodynamics/Hydrodynamics.cc`](https://github.com/EthanNT1000/asv_wave_sim/blob/ament_environment_hooks/gz-waves/src/systems/hydrodynamics/Hydrodynamics.cc)

**Primary Entry Point**: `Hydrodynamics::Update(...)` (`Physics.cc`)

**License**: GNU General Public License v3.0

## 1. Purpose and Scope

This implementation provides a physics-based model for computing the net hydrodynamic (and, optionally, aerodynamic) forces and torques acting on a rigid-body link (typically a marine vehicle hull represented as a triangular mesh) immersed in a dynamic wavefield and a bulk water current, within the Gazebo simulation framework. The model accounts for:

- **Buoyancy** (restoring force due to displaced fluid, computed at the center of pressure).
- **Viscous drag** (skin-friction drag derived from the 1957 ITTC correlation line).
- **Pressure drag** (form drag, differentiated for positive pressure and suction).
- **Per-DOF (Fossen-style) linear and quadratic damping** on surge/sway/heave/roll/pitch/yaw.
- **Planing-hull foil lift** on near-horizontal, bottom-facing submerged triangles.
- **Above-waterline aerodynamic drag** from ambient wind (flat-plate approximation).
- **A bulk water-current field**, sampled from a pre-processed HEC-RAS grid and added to the wave orbital velocity to give the total fluid velocity seen by every force term.

All geometric calculations operate in the CGAL kernel for exact triangle handling, vertex sorting, and vector operations. The wavefield interaction is mediated through a `WavefieldSampler` that supplies local water-surface heights, depths, and orbital velocities. The submerged-triangle and per-triangle force loops are parallelized with OpenMP.

## 2. High-Level Algorithm Flow

The algorithm executes once per simulation update via `Hydrodynamics::Update(...)`, which receives:

- Current wavefield sampler.
- Pose of the link's center of mass (CoM).
- Linear velocity $\mathbf{v}$ and angular velocity $\boldsymbol{\omega}$ of the CoM.
- The simulation time (`simTime`), used to sample the instantaneous wave orbital velocity.

The sequence is:

1. **Update submerged geometry** (`UpdateSubmergedTriangles`) — determine which triangles (or portions thereof) lie below the instantaneous water surface. Runs over mesh faces in parallel (OpenMP), with per-thread scratch buffers merged serially afterwards.
2. **Compute global quantities** — total hull area and submerged area (`ComputeAreas`), waterline length (`ComputeWaterlineLength`), waterline beam (`ComputeWaterlineBeam`), and a dynamic foil aspect ratio (`ComputeDynamicFoilGeometry`) used by the lift model.
3. **Sample the bulk current at the CoM** (`SampleWaterCurrentCoM`) — used by the Reynolds number and the per-DOF damping term.
4. **Single-pass per-triangle force loop** (`ComputeAllSubmergedForces`) — for every submerged triangle, in one parallel (OpenMP reduction) loop:
   - Point velocities and the relative-to-fluid velocity decomposition (`ComputePointVelocities`).
   - Buoyancy (center-of-pressure formulation).
   - Viscous drag (Reynolds-number dependent), if enabled.
   - Pressure drag (incidence-angle dependent), if enabled.
   - Foil lift (planing-hull dynamic lift), if enabled.
5. **Per-DOF damping** (`ComputeDampingForce`) — a single whole-body force/torque computed once from the CoM velocity relative to the bulk current, if enabled.

Forces and torques are accumulated at the CoM and returned for application to the rigid-body dynamics solver. A separate, optional above-waterline aerodynamic drag pass runs in `Hydrodynamics.cc` (the Gazebo system, not `Physics.cc`) using the same per-triangle geometry — see §4.9.

## 3. Core Data Structures

- **`HydrodynamicsParameters`** – Configurable coefficients (per-DOF damping, viscous, pressure-drag, foil lift) loaded from SDF or protobuf messages, plus a handle to the loaded `WaterCurrentGrid`. Supports optional randomization (`<randomize>`) for Monte-Carlo / domain-randomization studies.
- **`TriangleProperties`** – Per-mesh-triangle metadata: normal, area, submerged area, vertex height map (relative to water surface), and sorted high/mid/low vertices.
- **`SubmergedTriangleProperties`** – Per-submerged-triangle metadata: centroid $\mathbf{r}$, relative position $\mathbf{x}_r = \mathbf{r} - \mathbf{x}_{\text{CoM}}$, hull point velocity $\mathbf{v}_p$, the fluid velocity decomposition ($\mathbf{v}_{\text{orbital}}$, $\mathbf{v}_{\text{current}}$, $\mathbf{v}_{\text{fluid}}$), the hull-relative-to-fluid velocity $\mathbf{v}_{\text{rel}}$ and its normal/tangential split, angle of attack $\alpha$, and the direction of tangential flow.
- **`WaterCurrentGrid`** (`gz-waves/include/gz/waves/WaterCurrentGrid.hh`) – Loads a pre-processed binary grid (`WCRG v1`) and bilinearly samples a 2D bulk current velocity at any `(x, y)`.
- **`HydrodynamicsPrivate`** – Internal state holding the mesh, sampler, current pose/velocities, waterline segments, dynamic foil aspect ratio, sampled water current at the CoM, accumulated force/torque, and the persistent per-thread scratch buffers used by the parallel loops.

## 4. Detailed Computational Steps

### 4.1 Submerged-Triangle Identification (`UpdateSubmergedTriangles`)

For each vertex, depth is computed once via `WavefieldSampler::ComputeDepth` into a cached per-vertex property map (added once in the constructor), then for each face of the link mesh:

- Read the height map (negative depth) for the face's three vertices.
- Sort vertices by height: $H$ (highest), $M$ (middle), $L$ (lowest).

Both the vertex-depth loop and the face loop are parallelized with `#pragma omp parallel for`, with the thread count auto-scaled to the mesh size (at least 32 faces/vertices per thread) so small meshes don't pay OpenMP's barrier overhead. Each thread writes into its own persistent scratch buffers (`tl_subTris`, `tl_subProps`, `tl_waterlines`), which are merged into the shared `submergedTriangles` / `submergedTriangleProperties` / `waterline` vectors in a serial step afterwards.

Three cases are handled per face:

- **Fully submerged** ($h_H \leq 0$): The original triangle is added unchanged.
- **Partially submerged (two vertices above water)**: Split into one submerged triangle using linear interpolation parameters
  $$
  t_m = \frac{-h_L}{h_M - h_L}, \quad t_h = \frac{-h_L}{h_H - h_L}.
  $$
- **Partially submerged (one vertex above water)**: Split into two submerged triangles using analogous interpolation.

Orientation is preserved by ensuring the triangle normal matches the original face normal. Waterline segments are recorded for later length/beam computation, and each face's `subArea` is written back to `TriangleProperties` — the above-waterline aerodynamic drag pass (§4.9) uses `area - subArea` to get the exposed area of partially-submerged faces.

### 4.2 Global Geometric Quantities

- **Areas** (`ComputeAreas`):
  $$
  A_{\text{total}} = \sum \text{area}_i, \quad A_{\text{submerged}} = \sum \text{area}_j^{\text{sub}}.
  $$
- **Waterline length** (`ComputeWaterlineLength`): Projection of waterline edges onto the body-fixed x-axis, halved to avoid double-counting.
  $$
  L_{\text{waterline}} = \frac{1}{2} \sum_k \left| \mathbf{l}_k \cdot \hat{\mathbf{x}} \right|.
  $$
- **Waterline beam** (`ComputeWaterlineBeam`): the lateral (body-fixed y-axis) extent of all waterline segment endpoints,
  $$
  B_{\text{waterline}} = \max_k(\mathbf{l}_k \cdot \hat{\mathbf{y}}) - \min_k(\mathbf{l}_k \cdot \hat{\mathbf{y}}).
  $$
- **Dynamic foil aspect ratio** (`ComputeDynamicFoilGeometry`): "bottom" triangles are those with a submerged, upward-facing outward normal ($n_z >$ `BOTTOM_THRESHOLD`, i.e. within 72° of vertical). Their summed area gives the wetted bottom area $A_{\text{bottom}}$, and the aspect ratio used by the foil lift model is
  $$
  AR = \frac{B_{\text{waterline}}^2}{A_{\text{bottom}} + \epsilon}.
  $$

### 4.3 Point Velocities and Relative Fluid Velocity (`ComputePointVelocities`)

For each submerged triangle, given the CoM position/velocities, the wavefield sampler, the simulation time `t`, and the water current grid:
$$
\mathbf{x}_r = \mathbf{r} - \mathbf{x}_{\text{CoM}}, \quad \mathbf{v}_p = \mathbf{v} + \boldsymbol{\omega} \times \mathbf{x}_r.
$$
The fluid velocity at the centroid is the sum of the instantaneous wave orbital velocity and the bulk current, both sampled at the centroid:
$$
\mathbf{v}_{\text{orbital}} = \text{WavefieldSampler::ComputeOrbitalVelocity}(\mathbf{r}, t), \quad
\mathbf{v}_{\text{current}} = \text{WaterCurrentGrid::SampleAt}(r_x, r_y),
$$
$$
\mathbf{v}_{\text{fluid}} = \mathbf{v}_{\text{orbital}} + \mathbf{v}_{\text{current}}, \quad
\mathbf{v}_{\text{rel}} = \mathbf{v}_p - \mathbf{v}_{\text{fluid}}.
$$
$\mathbf{v}_{\text{rel}}$ — the hull velocity *relative to the fluid* — replaces the raw point velocity $\mathbf{v}_p$ as the input to every downstream force term (viscous drag, pressure drag, foil lift). It is decomposed into components normal and tangential to the triangle:
$$
(\mathbf{v}_{\text{rel}} \cdot \mathbf{n})\,\mathbf{n} = \mathbf{v}_{\text{rel},n}, \quad \mathbf{v}_{\text{rel},t} = \mathbf{v}_{\text{rel}} - \mathbf{v}_{\text{rel},n},
$$
and the angle of attack used by the foil lift model is
$$
\alpha = \operatorname{atan2}\!\left(\mathbf{v}_{\text{rel}} \cdot \mathbf{n},\; |\mathbf{v}_{\text{rel},t}|\right).
$$
As before, a unit tangential-flow direction $\mathbf{u}_f = -\mathbf{v}_{\text{rel},t}/|\mathbf{v}_{\text{rel},t}|$ and $\mathbf{v}_f = |\mathbf{v}_{\text{rel}}|\,\mathbf{u}_f$ are also derived, for use by the viscous drag term.

The Reynolds number (`ComputeReynoldsNumber`) and the per-DOF damping term both use the *whole-body* relative velocity $\mathbf{v} - \mathbf{v}_{\text{current,CoM}}$, where $\mathbf{v}_{\text{current,CoM}}$ is the current sampled once at the CoM (`SampleWaterCurrentCoM`) rather than per-triangle.

### 4.4 Buoyancy Force (`ComputeBuoyancyForce`)

Buoyancy is computed per submerged triangle using the **center-of-pressure** formulation (`Physics::BuoyancyForceAtCenterOfPressure`), unchanged from the original model, and returns a `(force, torque)` pair that the caller accumulates directly (rather than mutating shared state):

1. Compute centroid depth $h_C = $ `ComputeDepth(C)`.
2. Sort vertices $H, M, L$ by height.
3. Split the triangle horizontally at the waterline into upper ($HMD$) and lower ($LMD$) sub-triangles, where $D$ is the horizontal intercept and $B$ is the midpoint between $M$ and $D$.
4. For each sub-triangle, compute the hydrostatic pressure force at its centroid and the center of pressure using analytic formulas for triangular pressure distributions (apex-up or apex-down).
5. Combine into an overall center of application via weighted interpolation (`CenterOfForce`), giving the resultant force $\mathbf{f} = \mathbf{n}(f_U + f_L)$ and torque $\boldsymbol{\tau} = \mathbf{x}_r \times \mathbf{f}$.

(See §4 of the previous revision, or the annotated code in [`Hydrodynamics physics functions.md`](Hydrodynamics%20physics%20functions.md), for the full derivation — this part of the algorithm is unchanged.)

### 4.5 Viscous Drag Force (`ComputeViscousDragForce`)

Skin-friction drag follows the flat-plate (1957 ITTC) approximation:
$$
C_F = \frac{0.075}{(\log_{10} R_n - 2)^2}, \quad R_n = \frac{|\mathbf{v} - \mathbf{v}_{\text{current,CoM}}|\, L}{\nu},
$$
where $L$ is waterline length and $\nu$ is kinematic viscosity ($R_n \geq 10^3$ enforced).

Per submerged triangle, using the relative-to-fluid speed $|\mathbf{v}_{\text{rel}}|$ in place of the old $|\mathbf{v}_p|$:
$$
f_{\text{drag}} = \frac{1}{2} \rho\, C_F\, A\, |\mathbf{v}_{\text{rel}}|, \qquad \mathbf{f} = \mathbf{v}_f\, f_{\text{drag}}, \qquad \boldsymbol{\tau} = \mathbf{x}_r \times \mathbf{f}.
$$

### 4.6 Pressure Drag Force (`ComputePressureDragForce`)

Normal-force model differentiated by incidence angle, using $v = |\mathbf{v}_{\text{rel}}| / v_R$ (reference speed scaling):

- **Positive pressure** ($\cos\theta \geq 0$): $F = -(c_{P1} v + c_{P2} v^2)\, S\, (\cos\theta)^{f_P}$.
- **Suction** ($\cos\theta < 0$): $F = (c_{S1} v + c_{S2} v^2)\, S\, (-\cos\theta)^{f_S}$.

$S$ is triangle area and coefficients $c_{P/S}$ and exponents $f_{P/S}$ are user-configurable. Force acts along the triangle normal; torque follows from the lever arm $\mathbf{x}_r$.

### 4.7 Per-DOF (Fossen) Damping (`ComputeDampingForce`)

Unlike the earlier scalar linear/angular model, damping is now applied **per degree of freedom**, in the body frame, following a Fossen-style rigid-body damping matrix. The CoM linear/angular velocity (relative to the bulk current for the linear term) is rotated into the body frame:
$$
\mathbf{v}_B = R^{-1}\bigl(\mathbf{v} - \mathbf{v}_{\text{current,CoM}}\bigr), \qquad \boldsymbol{\omega}_B = R^{-1}\boldsymbol{\omega},
$$
where $R$ is the world-to-body rotation of the CoM pose. Each of the six DOFs (surge $u$, sway $v$, heave $w$, roll $p$, pitch $q$, yaw $r$) gets an independent linear + quadratic damping coefficient pair, e.g. for surge:
$$
F_u = -\bigl(c_{U1}\, u + c_{U2}\, |u|\, u\bigr),
$$
and analogously $F_v, F_w$ (force) and $\tau_p, \tau_q, \tau_r$ (torque) using $c_{V*}, c_{W*}, c_{P*}, c_{Q*}, c_{N*}$. The resulting body-frame force/torque vectors are rotated back into the world frame before being accumulated. Unlike the buoyancy/drag terms, this is computed **once per link per step** (not per submerged triangle) and is **not** scaled by the submerged-area ratio — the coefficients themselves are expected to characterize the whole hull.

### 4.8 Planing-Hull Foil Lift (`ComputeFoilLiftForce`)

An optional dynamic-lift term (`<foil_lift_on>`) models the hull's bottom surface as a low-aspect-ratio foil, active only on "bottom" triangles ($n_z >$ `BOTTOM_THRESHOLD`) with non-negligible relative flow ($|\mathbf{v}_{\text{rel}}| > 10^{-4}$ m/s):

$$
C_l = \begin{cases}
C_{l,\alpha}\, \alpha & |\alpha| < \alpha_{\text{stall}} \\
C_{l,\max}\, \operatorname{sign}(\alpha) & \text{otherwise}
\end{cases}, \qquad C_{l,\alpha} = 2\pi\, c_{\text{Lift1}},
$$
$$
C_{di} = \frac{C_l^2}{\pi\, AR + \epsilon}, \qquad q_A = \frac{1}{2}\rho\, |\mathbf{v}_{\text{rel}}|^2\, S,
$$

where $AR$ is the dynamic foil aspect ratio from §4.2 (floored at 0.5). The lift direction is the component of the triangle normal perpendicular to the flow direction $\hat{\mathbf{u}}_p = \mathbf{v}_{\text{rel}}/|\mathbf{v}_{\text{rel}}|$:
$$
\hat{\mathbf{l}} = \frac{\mathbf{n} - (\mathbf{n}\cdot\hat{\mathbf{u}}_p)\,\hat{\mathbf{u}}_p}{|\mathbf{n} - (\mathbf{n}\cdot\hat{\mathbf{u}}_p)\,\hat{\mathbf{u}}_p|}, \qquad
\mathbf{F}_{\text{foil}} = \hat{\mathbf{l}}\,(C_l\, q_A) - \hat{\mathbf{u}}_p\,(C_{di}\, q_A),
$$
with torque $\boldsymbol{\tau} = \mathbf{x}_r \times \mathbf{F}_{\text{foil}}$. The stall cap on $C_l$ prevents runaway lift at large angles of attack, and the induced-drag term $C_{di}$ approximates the trailing-vortex penalty of a finite-aspect-ratio foil.

### 4.9 Above-Waterline Aerodynamic Drag (`Hydrodynamics.cc`, not `Physics.cc`)

When `<aerodynamic_drag_on>` is set, a flat-plate wind-pressure model runs over each hull's `TriangleProperties` (already computed by `UpdateSubmergedTriangles`) in a separate OpenMP-parallel reduction loop in the Gazebo system, using the world `Wind` entity's linear velocity as $\mathbf{v}_{\text{wind}}$:
$$
A_{\text{above}} = \begin{cases} \text{area} & \text{fully above water} \\ \text{area} - \text{subArea} & \text{partially submerged} \end{cases},
$$
$$
v_n = (\mathbf{v}_{\text{wind}} - \mathbf{v}_{\text{hull}}) \cdot \hat{\mathbf{n}}, \qquad
\mathbf{F} = \tfrac{1}{2}\,\rho_{\text{air}}\, c_{\text{AeroDrag}}\, A_{\text{above}}\, v_n^2\, \hat{\mathbf{n}} \quad (v_n > 0),
$$
with $\rho_{\text{air}} = 1.225\ \text{kg/m}^3$ and faces on the lee side ($v_n \leq 0$) skipped. Torque is $\boldsymbol{\tau} = \mathbf{r} \times \mathbf{F}$ about the CoM. The resulting force/torque are applied to the link in addition to (not instead of) the hydrodynamic force/torque from `Hydrodynamics::Update`.

### 4.10 Bulk Water Current (`WaterCurrentGrid`)

`<water_current_grid>` (inside `<hydrodynamics>`) optionally points to a binary grid file (format `WCRG v1`: a small header — magic, version, `nx`, `ny`, grid origin and cell size — followed by interleaved `(vx, vy)` `float` pairs, row-major), typically produced from HEC-RAS output by an external `preprocess_hecras.py` script (not part of this repository). `WaterCurrentGrid::SampleAt(x, y)` bilinearly interpolates the 2D current velocity at any world `(x, y)`, returning zero outside the grid extent (no extrapolation). This same current is what §4.3 adds to the wave orbital velocity, and what §4.5/4.7 subtract from the CoM velocity for the Reynolds number and per-DOF damping.

## 5. Parameterization and Extensibility

All drag, damping, lift and current-related coefficients are exposed through `HydrodynamicsParameters::SetFromSDF` and `SetFromMsg`:

- Per-DOF damping: `cDampU1/U2`, `cDampV1/V2`, `cDampW1/W2`, `cDampP1/P2`, `cDampQ1/Q2`, `cDampN1/N2`.
- Pressure drag: `cPDrag1/2`, `fPDrag`, `cSDrag1/2`, `fSDrag`, `vRDrag` (unchanged).
- Foil lift: `foil_lift_on`, `cLift1`, `alphaStall`, `cLMax` (defaults to `cLift1 · 2π · sin(alphaStall)` if not set explicitly).
- Above-waterline drag: `aerodynamic_drag_on`, `cAeroDrag` (top-level plugin elements, not inside `<hydrodynamics>`).
- Water current: `water_current_grid` (path to a `WCRG v1` binary file).

An optional `<randomize>` SDF block draws the per-DOF damping and pressure-drag coefficients from uniform distributions (seeded from the system clock), facilitating sensitivity analysis and domain randomization; see `HydrodynamicsParameters::SetRandomFromSDF`.

Physical constants ($\rho$, $g$, $\nu$) reside in `PhysicalConstants.cc`:

- Gravity: -9.81, uniform acceleration due to gravity at earth's surface (orientation is z-up).
- G: 6.67408E-11 [m³ kg⁻¹ s⁻²], universal gravitational constant.
- WaterDensity: 1025 kg/m³.
- WaterKinematicViscosity: 1.0533E-6 m²/s.
- Air density for aerodynamic drag (`Hydrodynamics.cc`, not `PhysicalConstants.cc`): 1.225 kg/m³.

## 6. Implementation Notes and Optimizations

- **Parallelism**: `UpdateSubmergedTriangles` parallelizes the per-vertex depth computation and the per-face submersion/splitting loop; `ComputeAllSubmergedForces` parallelizes the per-submerged-triangle force loop (point velocities + buoyancy + viscous drag + pressure drag + foil lift) with an OpenMP reduction over the six scalar force/torque components. Thread counts are capped so each thread gets a minimum amount of work (≥32 faces/triangles), and per-thread scratch buffers are sized once and reused across steps to avoid per-tick allocation.
- **Numerical robustness**: Division-by-zero safeguards (e.g. `CenterOfForce`, foil aspect ratio, induced drag) and Reynolds-number clamping are present.
- **Debugging**: `DebugPrint` helpers and (mostly commented-out) `gzmsg` blocks allow inspection of intermediate triangle properties, forces, and centers of pressure.
- **Telemetry**: per-triangle and per-submerged-triangle properties can be streamed as InfluxDB line-protocol over UDP (`<influxDBUdp>` on the Gazebo plugin), and the model's speed through water (hull velocity relative to wave + current) is published on a Gazebo topic (`<SpeedThroughWater>`).
- **Limitations**: The model assumes a rigid, non-deforming mesh and neglects wave-making resistance or free-surface elevation feedback to the wavefield solver.
