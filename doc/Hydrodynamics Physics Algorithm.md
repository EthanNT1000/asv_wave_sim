# Hydrodynamics Physics Algorithm

**Source File**: [`gz-waves/src/Physics.cc`](https://github.com/EthanNT1000/asv_wave_sim/blob/ament_environment_hooks/gz-waves/src/Physics.cc)

**Primary Entry Point**: `Hydrodynamics::Update(...)` (line 811)

**License**: GNU General Public License v3.0

## 1. Purpose and Scope

This implementation provides a physics-based model for computing the net hydrodynamic forces and torques acting on a rigid-body link (typically a marine vehicle hull represented as a triangular mesh) immersed in a dynamic wavefield within the Gazebo simulation framework. The model accounts for:

- **Buoyancy** (restoring force due to displaced fluid, computed at the center of pressure).
- **Viscous drag** (skin-friction drag derived from the 1957 ITTC correlation line).
- **Pressure drag** (form drag, differentiated for positive pressure and suction).
- **Linear and quadratic damping** (scaled by the submerged-area ratio).

All calculations operate in the CGAL geometric kernel for exact triangle handling, vertex sorting, and vector operations. The wavefield interaction is mediated through a `WavefieldSampler` that supplies local water-surface heights and depths.

## 2. High-Level Algorithm Flow

The algorithm executes once per simulation update via `Hydrodynamics::Update(...)`, which receives:

- Current wavefield sampler.
- Pose of the link’s center of mass (CoM).
- Linear velocity $\mathbf{v}$ and angular velocity $\boldsymbol{\omega}$ of the CoM.

The sequence is:

1. **Update submerged geometry**
   Determine which triangles (or portions thereof) lie below the instantaneous water surface.
2. **Compute global quantities**
   Total hull area, submerged area, and waterline length.
3. **Compute point velocities**
   Velocity at each submerged-triangle centroid.
4. **Compute forces and torques**
   - Buoyancy (center-of-pressure formulation).
   - Viscous drag (Reynolds-number dependent).
   - Pressure drag (incidence-angle dependent).
   - Damping (linear + quadratic, area-ratio scaled).

Forces and torques are accumulated at the CoM and returned for application to the rigid-body dynamics solver.

## 3. Core Data Structures

- **`HydrodynamicsParameters`** – Configurable coefficients (damping, viscous, pressure-drag) loaded from SDF or protobuf messages. Supports optional randomization for Monte-Carlo studies.
- **`TriangleProperties`** – Per-mesh-triangle metadata: normal, area, vertex height map (relative to water surface), and sorted high/mid/low vertices.
- **`SubmergedTriangleProperties`** – Per-submerged-triangle metadata: centroid $\mathbf{r}$, relative position $\mathbf{x}_r = \mathbf{r} - \mathbf{x}_{\text{CoM}}$, velocity decomposition into normal/tangential components, and direction of tangential flow.
- **`HydrodynamicsPrivate`** – Internal state holding the mesh, sampler, current pose/velocities, waterline segments, and accumulated force/torque.

## 4. Detailed Computational Steps

### 4.1 Submerged-Triangle Identification (`UpdateSubmergedTriangles`)

For each face of the link mesh:

- Compute the signed depth of each vertex using `WavefieldSampler::ComputeDepth`.
- Store the height map (negative depth) in `TriangleProperties::heightMap`.
- Sort vertices by height: $H$ (highest), $M$ (middle), $L$ (lowest).

Three cases are handled:

- **Fully submerged** ($h_H \leq 0$): The original triangle is added unchanged.
- **Partially submerged (two vertices above water)**: Split into one submerged triangle using linear interpolation parameters
  $$
  t_m = \frac{-h_L}{h_M - h_L}, \quad t_h = \frac{-h_L}{h_H - h_L}.
  $$
- **Partially submerged (one vertex above water)**: Split into two submerged triangles using analogous interpolation.

Orientation is preserved by ensuring the triangle normal matches the original face normal. Waterline segments are recorded for later length computation.

Key Definitions:

- Height map: $  \text{heightMap}[i] = - \text{depth}(\text{vertex}_i)  $
- Sorted indices: $  \text{idx} = \text{algorithm::sort\_indexes}(\text{heightMap})  $
- Submerged area fraction: $  \text{subArea} = \text{area of generated submerged triangle(s)}  $

### 4.2 Global Geometric Quantities

- **Areas** (`ComputeAreas`):
  $$
  A_{\text{total}} = \sum \text{area}_i, \quad A_{\text{submerged}} = \sum \text{area}_j^{\text{sub}}.
  $$
- **Waterline length** (`ComputeWaterlineLength`): Projection of waterline edges onto the body-fixed x-axis, halved to avoid double-counting.
  $$
  L_{\text{waterline}} = \frac{1}{2} \sum_k \left| \mathbf{l}_k \cdot \hat{\mathbf{x}} \right|,
  $$
  where $  \mathbf{l}_k  $ is the vector of each waterline segment and $  \hat{\mathbf{x}}  $ is the body-fixed x-direction.

### 4.3 Point Velocities at Centroids (`ComputePointVelocities`)

For each submerged triangle:
$$
\mathbf{x}_r = \mathbf{r} - \mathbf{x}_{\text{CoM}}, \quad \mathbf{v}_p = \mathbf{v} + \boldsymbol{\omega} \times \mathbf{x}_r.
$$
Decompose:
$$
\mathbf{u}_p = \frac{\mathbf{v}_p}{|\mathbf{v}_p|}, \quad \cos\theta = \mathbf{u}_p \cdot \mathbf{n},
$$
$$
\mathbf{v}_n = (\cos\theta)\,\mathbf{n}, \quad \mathbf{v}_t = \mathbf{v}_p - \mathbf{v}_n,
$$
$$
\mathbf{u}_f = -\frac{\mathbf{v}_t}{|\mathbf{v}_t|}, \quad \mathbf{v}_f = |\mathbf{v}_p|\,\mathbf{u}_f.
$$

Key Definitions:

- $ \mathbf{x}_r $: relative position of the centroid wrt CoM
- $  \mathbf{r}  $: triangle centroid
- $  \mathbf{v} $: linear velocity of the centre of mass
- $  \boldsymbol{\omega} $: angular velocity of the centre of mass
- $  \mathbf{v}_p  $: velocity at centroid (translational + rotational contribution)
- $  \mathbf{u}_p  $: normalized point velocity
- $  \mathbf{n}  $: triangle normal
- $  \cos\theta  $: cosine of angle between velocity and surface normal
- $  \mathbf{v}_n  $: point velocity normal to surface
- $  \mathbf{v}_t  $: point velocity tangential to surface
- $  \mathbf{u}_f  $: direction of tangential flow
- $  \mathbf{v}_f  $: tangential flow velocity vector used for viscous drag

### 4.4 Buoyancy Force (`ComputeBuoyancyForce`)

Buoyancy is computed per submerged triangle using the **center-of-pressure** formulation (`Physics::BuoyancyForceAtCenterOfPressure`) rather than the simple centroid method, providing higher fidelity for sloped or partially submerged faces.

**Inner algorithm** (`BuoyancyForceAtCenterOfPressure` overload):

1. Compute centroid depth $h_C = $ `ComputeDepth(C)`.
1. Sort vertices $H, M, L$ by height.
1. Split the triangle horizontally at the waterline to obtain upper ($HMD$) and lower ($LMD$) sub-triangles, where $D$ is the horizontal intercept and $B$ is the midpoint between $M$ and $D$.
1. For each sub-triangle, calculate:
   - Hydrostatic pressure force at its centroid:
     $$
      f_U = \rho g A_{HMD} h_{CU}, \quad f_L = \rho g A_{LMD} h_{CL}
     $$
     where $\rho$ is water density and $g$ is gravity.
   - Center of pressure using analytic formulas for triangular pressure distributions (apex-up or apex-down):
      $$
        t_{cU} = \frac{4z_0 + 3h}{6z_0 + 4h}
      $$
      $$
        \mathbf{C}_{pU} = H + alt(B - H) \times t_{cU} \quad (\text{if } H.z > M.z)
      $$
      $$
        t_{cL} = \frac{2z_0 + h}{6z_0 + 2h}
      $$
      $$
        \mathbf{C}_{pL} = B + alt(L - B) \times _{cL} \quad (\text{if } M.z > L.z)
      $$

      - $  \mathbf{C}_{pU}, \mathbf{C}_{pL}  $: centers of pressure for upper/lower sub-triangles
      - $  z_0  $: reference depth offset for center-of-pressure calculation
1. Combine forces and locate the overall center of application via weighted interpolation (`CenterOfForce`):

  $$
  \mathbf{C_f} = \frac{f_U \mathbf{C}_{pU} + f_L \mathbf{C}_{pL}}{f_U + f_L} \quad (\text{if } |f_U + f_L| > \epsilon)
  $$
  $$
  \mathbf{f} = n (f_U + f_L)
  $$
  $$
  \boldsymbol{\tau} = \mathbf{x}_r \times \mathbf{f}
  $$

- $\mathbf{C}$ : center of force
- $\mathbf{f}$ : total buoyancy force vector
- $\mathbf{n}$ : normalized triangle
- $\boldsymbol{\tau}$ : torque

### 4.5 Viscous Drag Force (`ComputeViscousDragForce`)

Skin-friction drag follows the flat-plate approximation:
$$
C_F = \frac{0.075}{(\log_{10} R_n - 2)^2}, \quad R_n = \frac{u L}{\nu},
$$
where $L$ is waterline length, $u = |\mathbf{v}|$, and $\nu$ is kinematic viscosity (with $R_n \geq 10^3$ enforced).

Per submerged triangle:
$$
f_{\text{drag}} = \frac{1}{2} \rho C_F A |\mathbf{v}_f|,
$$
Applied in the direction of tangential flow $\mathbf{v}_f$:
$$
force = vf \times f_{drag}
$$
$$
\boldsymbol{\tau} = \mathbf{x}_r \times force
$$

### 4.6 Pressure Drag Force (`ComputePressureDragForce`)

Normal-force model differentiated by incidence angle:

- **Positive pressure** ($\cos\theta \geq 0$):
  $$
  F = -(c_{P1} v + c_{P2} v^2) S (\cos\theta)^{f_P}.
  $$
- **Suction** ($\cos\theta < 0$):
  $$
  F = (c_{S1} v + c_{S2} v^2) S (-\cos\theta)^{f_S}.
  $$

Here $v = v_p / v_{R}$ (reference speed scaling), $S$ is triangle area, and coefficients $c_{P/S}$ and exponents $f_{P/S}$ are user-configurable.

Force acts along the triangle normal; torque follows from the lever arm $\mathbf{x}_r$.

### 4.7 Damping Force (`ComputeDampingForce`)

A simple area-ratio-scaled model:
$$
r_s = \frac{A_{\text{submerged}}}{A_{\text{total}}},
$$
$$
\mathbf{f} = r_s \bigl( -c_{L1} - c_{L2} |\mathbf{v}| \bigr) \mathbf{v},
$$
$$
\boldsymbol{\tau} = r_s \bigl( -c_{R1} - c_{R2} |\boldsymbol{\omega}| \bigr) \boldsymbol{\omega}.
$$
This term provides numerical stability and low-speed dissipation.

## 5. Parameterization and Extensibility

All drag and damping coefficients are exposed through `HydrodynamicsParameters::SetFromSDF` and `SetFromMsg`. An optional `<randomize>` SDF block draws coefficients from uniform distributions, facilitating sensitivity analysis.

Physical constants ($\rho$, $g$, $\nu$) reside in `PhysicalConstants.cc`.

- Gravity: -9.81, Uniform acceleration due to gravity at earth's surface (orientation is z-up)
- G : 6.67408E-11 [m3 kg-1 s-2], Universal gravitational constant.
- WaterDensity : 1025
- WaterKinematicViscosity: 1.0533E-6

## 6. Implementation Notes and Optimizations

- **Numerical robustness**: Division-by-zero safeguards (e.g., `CenterOfForce`) and Reynolds-number clamping are present.
- **Performance**: Vertex depths are cached once per update; triangle splitting reuses geometry utilities from `Geometry.hh`.
- **Debugging**: Extensive commented `gzmsg` blocks and `DebugPrint` helpers allow inspection of intermediate triangle properties, forces, and centers of pressure.
- **Limitations**: The model assumes a rigid, non-deforming mesh and neglects wave-making resistance or free-surface elevation feedback to the wavefield solver.
