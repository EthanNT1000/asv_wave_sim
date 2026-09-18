# Hydrodynamics audit against Fossen (Handbook of Marine Craft Hydrodynamics and Motion Control, 2nd ed.)

Scope: `gz-waves/src/Physics.cc`, `gz-waves/include/gz/waves/Physics.hh`, the
Gazebo system that calls them (`gz-waves/src/systems/hydrodynamics/Hydrodynamics.cc`),
the wave kinematics feeding them (`WavefieldSampler.cc`, `WaveParameters.cc`,
`LinearRegularWaveSimulation.cc`, `LinearRandomFFTWaveSimulation.cc`,
`WaveSpectrum.cc`, `WaveSpreadingFunction.cc`) and the current field
(`WaterCurrentGrid.cc`).

Reference docs available in the repo: `doc/Hydrodynamics Physics Algorithm.md`
and `doc/Hydrodynamics physics functions.md`. The two Word documents named in
the task (`wave_foil_dynamics_equations_v2.docx`,
`hydrodynamics_parameter_guide_v3.docx`) are **not present in the repository**;
this audit is therefore against Fossen and the code itself.

Citation convention: Fossen is cited by chapter/section and topic. Equation
numbers are given where they are stable across the two editions (e.g. Theorem
3.2, Eq. 2.18, Eq. 3.44); otherwise the section is the authoritative pointer.

Numerical checks quoted below were made with the pure-numpy port of the force
model in `tests/hydro/hydro_reference.py`, which reproduces `Physics.cc`
term by term (see `tests/hydro/test_properties.py`).

---

## 0. Executive summary

| ID | Term | Severity | Status |
|----|------|----------|--------|
| F-01 | No added mass `M_A`, no `C_A(nu_r)`, no Munk moment; hybrid model class | critical | not applied - needs decision |
| F-02 | Foil lift acts along `+n` (should be `-n`) and its "bottom face" gate selects upward-facing faces | critical | not applied - needs decision (patch provided) |
| F-03 | Orbital velocity is computed from the regular/trochoid component parameters even when the wavefield is the FFT sea (default) | critical | not applied - needs decision |
| F-04 | Wave loads are hydrostatic Froude-Krylov only: `p = rho g (eta - z)` with no `e^{kz}` decay, no diffraction, no radiation, no 2nd-order drift | moderate | not applied - needs decision |
| F-05 | Link-origin velocity used as CoM velocity in `v_p = v + omega x r` (also in the aero pass) | moderate | **applied** |
| F-06 | Aerodynamic drag loads the lee faces instead of the windward faces | moderate | **applied** |
| F-07 | Double counting of quadratic cross-flow/normal drag: pressure drag + per-DOF `cDampV2/N2` + foil normal force on the same faces | moderate | not applied - tuning decision |
| F-08 | Whole-body `D(nu_r)` is not scaled by wetted fraction (acts when the hull is airborne) | moderate | not applied - needs decision |
| F-09 | Skin friction uses the full relative speed, not the tangential component | moderate | not applied - needs decision |
| F-10 | Spectrum is ECKV (`U10`, wave age), no `Hs/Tp` (ITTC/JONSWAP/Torsethaugen) interface; spreading OK | moderate | not applied |
| F-11 | Numerics: explicit forces with one-step velocity lag; stability constraint documented | moderate | documented |
| F-12 | Foil `C_L_alpha = 2 pi` ignores finite aspect ratio; `AR` uses total waterline beam (catamaran) | cosmetic | not applied |
| F-13 | Per-DOF damping frame taken from the inertial (CoM) pose, not the link pose | cosmetic | not applied |
| F-14 | `cLMax` default `2 pi sin(alpha_stall)` is 1 % discontinuous with the linear branch | cosmetic | not applied |
| F-15 | Spatially varying current grid violates Fossen's constant irrotational current assumption (`dot nu_c != 0`) | cosmetic | documented |
| F-16 | Restoring `g(eta)`: exact nonlinear hydrostatics, `GM_T`, `GM_L`, `A_wp` verified correct | none | verified |
| F-17 | Encounter frequency: not needed explicitly, emerges from evaluating the field at the body position | none | verified |
| F-18 | Relative velocity enters damping/drag only, not rigid-body Coriolis; rotational DOFs excluded from `nu_c` | none | verified |

---

## 1. Model class (Fossen Ch. 5-6)

Fossen's unified equation of motion for a vessel in a current (Sec. 6.3 and
Sec. 7.1, "maneuvering model in the presence of ocean currents"):

```
M_RB nu_dot + C_RB(nu) nu + M_A nu_r_dot + C_A(nu_r) nu_r + D(nu_r) nu_r + g(eta) = tau + tau_wind + tau_wave
```

with `M = M_RB + M_A` and, for seakeeping (Sec. 5.4, Cummins equation),
`M_A` replaced by `A(inf)` plus the retardation convolution `int K(t - tau) nu_r(tau) dtau`.

How the code maps onto this:

| Fossen term | Where it lives | Notes |
|---|---|---|
| `M_RB nu_dot + C_RB(nu) nu` | Gazebo physics engine (DART) from the SDF `<inertial>` | Not in `Physics.cc`. Engine integrates with the absolute velocity `nu`, which is correct (Fossen Sec. 7.1: the current is irrotational, so `C_RB` is written in `nu`, not `nu_r`). |
| `M_A nu_r_dot` | **absent** | see F-01 |
| `C_A(nu_r) nu_r` (incl. Munk moment) | **absent** | see F-01 |
| Retardation `K(t)` (Cummins) | **absent** | pure zero-memory model |
| `D_L nu_r` | `ComputeDampingForce` (`cDamp*1`) | diagonal, body frame, `nu_r = nu - nu_c` for the 3 linear DOFs |
| `D_n(nu_r) nu_r` | `ComputeDampingForce` (`cDamp*2`) plus per-triangle `ComputePressureDragForce`, `ComputeViscousDragForce`, induced drag in `ComputeFoilLiftForce` | see F-07 for overlap |
| `g(eta)` | `ComputeBuoyancyForce` (pressure integration over the instantaneous wetted surface) | exact nonlinear hydrostatics, verified (F-16) |
| `tau_wave` (1st order) | same buoyancy integration, because the depth is measured to the *local wave surface* | hydrostatic Froude-Krylov only (F-04) |
| `tau_wave` (2nd order drift) | **absent** | F-04 |
| `tau_wind` | aero pass in `Hydrodynamics.cc` | face-based, see F-06 |
| `nu_c` | `WaterCurrentGrid::SampleAt` (2-D, depth independent) | enters drag/damping only (F-18) |

**Verdict.** The implementation is neither a Fossen maneuvering model (no
zero-frequency added mass) nor a seakeeping model (no radiation memory). It is
a hydrostatic + empirical viscous model in the Kerner/SimShip tradition, to
which a diagonal `D_L + D_n` has been bolted on. Because the rigid-body part is
integrated by the engine and every hydrodynamic term is an external force
evaluated from the previous state, the model is internally *consistent* as a
zero-memory model, but it is missing the inertial hydrodynamic terms entirely.

### F-01 (critical) No added mass, no `C_A(nu_r)`, no Munk moment

* **As implemented:** `M = M_RB`. `Hydrodynamics::Update` produces only
  velocity- and position-dependent forces. Nothing depends on `nu_dot`.
* **Fossen:** Sec. 6.3, `M_A = -[X_udot ...]` symmetric positive (semi-)definite for
  a body at rest in an ideal fluid; `C_A(nu_r)` obtained from `M_A` via the
  parametrization of Theorem 3.2 (Sec. 6.3, "hydrodynamic Coriolis and centripetal
  matrix"); the destabilising Munk moment in pitch/yaw is the
  `(Z_wdot - X_udot) u_r w_r` / `(Y_vdot - X_udot) u_r v_r` component of `C_A(nu_r) nu_r`.
* **Discrepancy:** `M_A = 0`, so `C_A = 0` and the Munk moment is absent. For a
  barge-like or twin-pontoon hull the heave/pitch/roll added mass is of the same
  order as `M_RB`, so the natural periods are too short by roughly
  `sqrt((M_RB + M_A) / M_RB)` (20-40 %), accelerations in waves are overestimated, and
  the yaw/pitch behaviour at speed lacks the Munk destabilisation that the real
  hull has (the `cDampN*` values then have to be tuned to compensate).
* **Consequence for the audit questions:** `M` is trivially symmetric positive
  definite (it is the SDF inertia), and there is no hand-written `C_A` to
  check against a derived one; the skew-symmetry test in
  `tests/hydro/test_properties.py` therefore exercises the *reference*
  construction `C(nu) = f(M)` that any future `M_A` must go through.
* **Fix options (needs decision):**
  1. *Preferred, zero code:* declare `<fluid_added_mass>` inside `<inertial>` in
     the model SDF (SDF 1.10 / Gazebo Garden+, DART). The engine then adds
     `M_A` to the spatial inertia and its Coriolis term follows from the full
     inertia, i.e. `C_A` is derived, not hand-written, and the Munk moment
     appears automatically. Caveat: the engine uses `nu`, not `nu_r`, in
     `C_A`; the difference `C_A(nu) nu - C_A(nu_r) nu_r` is a current-dependent
     bias that is small for `|nu_c| << |nu|`.
  2. Apply `-M_A nu_dot - C_A(nu_r) nu_r` as an external force in the plugin.
     Not recommended: `nu_dot` must be estimated by finite differences from
     lagged velocities, which creates an algebraic loop that is unstable for
     `M_A` comparable to `M_RB` (see Sec. 7 below).
  3. Move to a Cummins formulation (Sec. 5.4) with `A(inf)` and a state-space
     approximation of `K(t)`. This is a modelling project (needs hydrodynamic
     coefficients from a panel code), not a patch.

### F-02 (critical) Foil lift direction and "bottom" gate are inconsistent with the mesh normal convention

* **As implemented** (`ComputeFoilLiftForce`, `ComputeDynamicFoilGeometry`):
  a face qualifies as "bottom" if its normal has `n_z > 0.3` (upward-facing);
  `alpha = atan2(v_rel . n, |v_rel_t|)` with `v_rel` the hull velocity relative to
  the fluid; `C_l = 2 pi cLift1 alpha`; force `F = lift_dir * C_l q A` with
  `lift_dir = normalize(n - (n . u) u)`, i.e. along `+n` for `alpha > 0`.
* **Normal convention required by the rest of the file:** the buoyancy term is
  `F = n * rho * g * A * h` with `g = -9.81`, i.e. `F = -rho |g| A h n`. This
  is only upward on the hull bottom if `n` is the *outward* normal. The
  shipped WAM-V collision meshes have positive signed volume (outward), and the
  pressure-drag term uses the same convention correctly (a face moving into the
  fluid, `v_rel . n > 0`, is pushed along `-n`).
* **Discrepancy 1 (gate):** with outward normals the planing bottom has
  `n_z < 0`. The `n_z > 0.3` test excludes the bottom and instead selects
  upward-facing faces (deck, pontoon tops) whenever a wave wets them.
* **Discrepancy 2 (sign):** for a bow-up planing bottom (outward normal down and
  slightly forward, hull moving forward) the port gives
  `alpha = +tau_trim > 0` and `F = (-36, 0, -725) N` for a 0.5 m^2 face at 3 m/s,
  i.e. **downward**; physically the flow impinges on the outside of the face and
  the resultant is along `-n` (upward). The same sign inversion appears on any
  face the gate does select. In general the normal force on a plate must oppose
  the normal relative velocity: `F ~ -sign(v_rel . n) n`, exactly as the
  pressure-drag term does; the foil term has the opposite sign.
* **Fossen reference:** lift on a foil/rudder, Ch. 9 (control surfaces),
  `L = 1/2 rho A C_L(alpha) V^2` directed perpendicular to the inflow toward
  the suction side; the sign convention follows the angle of attack measured
  from the inflow to the chord.
* **Why not applied:** the change flips the sign of an active default-on force
  term (`foilLiftOn(true)`) and changes which faces it acts on; the commit
  history (`3ab7ac7`) shows the current behaviour was tuned empirically
  ("sudden high thrust from rest now causes pitch-up flip"). This needs the hull
  owner to confirm their mesh orientation and re-tune. The fix, for an
  outward-normal mesh, is two lines:

  ```cpp
  // ComputeDynamicFoilGeometry and ComputeFoilLiftForce: bottom = downward-facing
  if (nz > -bottomThresh || props.v_rel_mag < 1e-4) return {NULL, NULL};   // was: nz < bottomThresh
  ...
  // resultant opposes the normal relative velocity (same convention as pressure drag)
  const cgal::Vector3 F_foil = lift_dir * (-Cl * q_A) + (-props.up) * (Cdi * q_A);
  ```
  (`ComputeDynamicFoilGeometry` must use the same `nz < -BOTTOM_THRESHOLD` test.)
  The test `test_foil_lift_sign_reference` in `tests/hydro/test_properties.py`
  encodes the expected behaviour and is marked expected-failure against the
  current formula.

---

## 2. `M = M_RB + M_A` and `C_A(nu)` derivation

* `M_RB` (Fossen Eq. 3.44) is built by the engine from `<mass>`, `<inertia>` and
  the inertial `<pose>`; the WAM-V model has `r_g = 0`, `m = 180 kg`,
  `I = diag(120, 393, 446)`. Symmetric positive definite by construction;
  asserted in `test_mass_matrix_spd`.
* `M_A`: absent (F-01). No hand-written `C_A` exists, so there is nothing to
  flag for "hand-written vs derived"; the reference test constructs
  `C(nu)` from `M` with Theorem 3.2 and asserts `C + C^T = 0`, including the
  Munk moment sign for a representative diagonal `M_A`.
* If `<fluid_added_mass>` is introduced (F-01 option 1) the SDF matrix must be
  symmetric and positive semi-definite; the test reads it from the SDF when
  present.

---

## 3. Damping `D(nu_r) = D_L + D_n(nu_r)` (Fossen Sec. 6.4)

* `ComputeDampingForce` implements exactly `tau_i = -(c1_i v_i + c2_i |v_i| v_i)`
  per DOF in the body frame, i.e. `D = diag(c1_i + c2_i |v_i|)`. This is the
  diagonal form of Fossen's `D_L + D_n(nu_r)` (Sec. 6.4.1 linear damping,
  Sec. 6.4.2 quadratic/cross-flow) and is positive definite whenever all
  `c1_i > 0` (asserted in `test_damping_positive_definite`). Energy:
  `nu_r^T D(nu_r) nu_r >= 0`, dissipative.
* Default and WAM-V values (`1e-6 ... 5e-3`) are physically negligible for a
  180 kg craft (linear surge damping of a WAM-V-class hull is O(10-100) N s/m);
  essentially all damping comes from the per-triangle terms. Not a bug, a
  tuning observation.

### F-07 (moderate) Double counting between foil model, pressure drag and per-DOF quadratic damping

Fossen Sec. 6.4.2 models quadratic sway/yaw damping as a *single* cross-flow
drag integral `Y = -1/2 rho int T(x) C_d^{2D}(x) |v_r + x r| (v_r + x r) dx` and
warns that when a strip/panel model provides it, the lumped `|v|v` terms must
not be added again. Here three terms provide a quadratic normal force on the
same wetted faces:

1. `ComputePressureDragForce`: `-(cPDrag1 v + cPDrag2 v^2) S cos^0.4(theta) n` is a
   distributed cross-flow drag (the `v^2 S cos theta` part is exactly
   `1/2 rho C_d S |v_n| v_n` with `C_d = 2 cPDrag2 / (rho vRDrag^2)`).
2. `cDampV2 |v| v`, `cDampN2 |r| r`: lumped cross-flow drag in sway and yaw.
3. `ComputeFoilLiftForce` on bottom faces: `C_l q A ~ 2 pi alpha 1/2 rho V^2 S`, with
   `alpha ~ v_n / V`, is again a force linear in `v_n` and in `V` on the same
   face that already carries the pressure-drag normal force.

Recommendation (tuning decision, not code): treat (1) as the cross-flow drag
and set `cDampV2 = cDampN2 = 0`, or keep (2) and set `cPDrag* = 0`; when the
foil term is enabled, reduce `cPDrag*` on the faces that pass the bottom gate.
The current defaults double count.

### F-08 (moderate) Whole-body damping is not scaled by wetted fraction

`ComputeDampingForce` applies the full `D(nu_r) nu_r` regardless of how much of
the hull is in the water, so a planing hull that leaves the water is still
damped in all 6 DOF. Fossen's `D` is defined for the floating vessel. Fix
(needs a reference wetted area): scale by `A_sub / A_sub,rest`, or gate on
`submergedArea > 0`.

### F-09 (moderate) Skin friction uses `|v_rel|`, not the tangential speed

`ComputeViscousDragForce`: `F = -u_t * 1/2 rho C_F A |v_rel|^2`. The ITTC-57
line (`C_F = 0.075 / (log10 Rn - 2)^2`, correctly implemented in
`ViscousDragCoefficient`) gives the shear from the *tangential* flow; the
correct magnitude is `1/2 rho C_F A |v_t|^2`. Using `|v_rel|^2` over-predicts
skin friction on faces whose motion is mostly normal (bottom in heave), where the
pressure-drag term already acts. This is Kerner's deliberate simplification
(documented in his article); changing it changes the heave/pitch damping and
requires re-tuning, so it is left as a decision.

---

## 4. Current / relative velocity (Fossen Sec. 7.1, Sec. 8.3)

Verified correct (F-18):

* `nu_c = [u_c v_c 0 0 0 0]` from `WaterCurrentGrid::SampleAt(x, y)` (2-D,
  `w_c = 0`). Rotational components are zero, matching Fossen's irrotational
  current assumption (`omega_c = 0`, so `omega_r = omega`).
* `nu_r` enters: per-triangle drag/lift (`v_rel = v_p - v_orbital - v_current`),
  the whole-body damping (`linVelocity - waterCurrentCoM`, linear DOFs only) and
  the Reynolds number. It does **not** enter buoyancy nor the rigid-body
  Coriolis term (which stays in the engine with `nu`). This is the energy
  consistent form: every term with `nu_r` is dissipative
  (`F . v_rel <= 0`, checked for each term in Sec. 6 of the reference port), and
  `C_RB(nu) nu` does no work.
* No `M_A`, so the `nu_r` vs `nu` question for added mass is moot (F-01).

F-15 (cosmetic): the HEC-RAS grid is spatially varying, so a moving hull sees
`d nu_c / dt != 0` and Fossen's "constant irrotational current in NED" model
(Sec. 8.3) is violated. The extra term `M_A nu_c_dot` would only matter once
`M_A` exists; the drag terms are unaffected.

### F-05 (moderate, applied) Link-origin velocity used as the CoM velocity

* **As implemented:** `Hydrodynamics.cc` passes `link.WorldLinearVelocity()`
  (velocity of the link origin) as the CoM velocity to `Hydrodynamics::Update`,
  which then forms `v_p = v + omega x (r - r_CoM)` (`ComputePointVelocities`)
  and `nu_r` for damping and Reynolds number. The aero pass does the same.
* **Fossen:** rigid-body kinematics, Sec. 3.1: `v_p = v_o + omega x r_p`
  measured from the *same* point whose velocity is used. Using the origin
  velocity with a CoM lever arm drops the term `omega x r_g`.
* **Effect:** zero when the inertial `<pose>` is at the link origin (WAM-V);
  otherwise every relative-velocity term (drag, lift, damping, Reynolds) is
  biased by `omega x r_g`, which mis-predicts roll/pitch damping torques.
* **Fix (applied):** `v_CoM = v_link + omega x (p_CoM - p_link)` computed once
  in the plugin and passed to `Update` and to the aero pass. No API change.

---

## 5. Wave loads (Fossen Sec. 5.1-5.3, Sec. 8.2)

### F-04 (moderate) Hydrostatic Froude-Krylov only

* **As implemented:** the only wave load is the buoyancy integral with the depth
  taken to the *local instantaneous* free surface: `p = rho g (eta(x,y,t) - z)`.
  Because it uses the exact wetted geometry it is a *nonlinear* Froude-Krylov
  force (better than Fossen's linearised `tau_wave1` for steep waves), but:
  * no `e^{kz}` attenuation and no dynamic pressure: the linear FK pressure is
    `p = rho g eta e^{kz} - rho g z` (Sec. 8.2, linear wave theory); for a
    WAM-V (draft ~0.2 m) in 2 s waves (`k ~ 1`) the code over-predicts the wave
    pressure on the bottom by ~20 %, and much more for short waves.
    `LinearRegularWaveSimulation::Pressure` already implements the correct
    `e^{kz}` scaling but is not used by `Physics.cc`.
  * no diffraction (valid only for `lambda >> L`; for a 4.9 m WAM-V that means
    `lambda > ~25 m`, i.e. `T > 4 s`), hence no force RAO;
  * no radiation (no `A(omega)`, `B(omega)`, see F-01);
  * no second-order (mean and slowly varying) drift forces (Sec. 8.2, wave
    drift), so station-keeping simulations have no wave drift load.
* **Fix (needs decision):** expose `Pressure()`/velocity from the wave
  simulations through `WavefieldSampler` and use `p = rho g (eta e^{kz} - z)` in
  `BuoyancyForceAtCenterOfPressure`; drift and diffraction require a separate
  force-RAO/QTF model (Sec. 5, 8.2), out of scope for a patch.

### F-03 (critical) Orbital velocity is inconsistent with the FFT wavefield

* `WavefieldSampler::ComputeOrbitalVelocity` sums linear deep-water orbital
  velocities (`u = a omega e^{kz} cos theta`, `w = a omega e^{kz} sin theta`,
  correct for `eta = a cos theta`) over `WaveParameters::Amplitude_V()` etc.
  Those component vectors are the regular/trochoid parameterisation
  (`number`, `scale`, `angle`, `steepness`).
* `WaveParameters::algorithm_` defaults to `"fft"`: the surface that lifts the
  hull is the ECKV random sea from `LinearRandomFFTWaveSimulation`, which does
  not use `Amplitude_V()` at all. The drag/lift terms therefore see orbital
  velocities of a *different* (regular) wave train than the one producing the
  buoyancy. For `linear_regular` and `trochoid` the two are consistent (phase
  convention `theta = k . x - omega t` verified in both).
* **Fix (needs decision):** add velocity spectra to the FFT simulation
  (`i omega h_k` terms, same machinery as the displacement FFTs) and route them
  through `WavefieldSampler`, or zero the orbital velocity when
  `algorithm == fft` so at least the model is consistent (still-water drag).

### Encounter frequency (verified, F-17)

`omega_e = omega - k U cos beta` (Fossen Sec. 8.2, encounter frequency) is not
applied explicitly and does not need to be: both the elevation and the orbital
velocity are evaluated at the body's actual `(x, y, t)`, so the Doppler shift
emerges from `theta = k . x(t) - omega t`. Correct.

### Spectrum and spreading (F-10, moderate)

* `LinearRandomFFTWaveSimulation` uses the ECKV (Elfouhaily et al. 1997)
  omnidirectional `k`-spectrum with `U10` and wave age `Omega_c`, and cos-2s
  spreading (`Cos2sSpreadingFunction`). `PiersonMoskowitzWaveSpectrum` exists
  and is a correct `k`-space transform of Fossen's PM `S(omega) = A omega^-5 exp(-B omega^-4)`
  (`alpha/(2k^3) exp(-beta g^2/(k^2 U^4))`, checked by substitution), but it is
  not wired to any algorithm.
* Fossen (Sec. 8.2, wave spectra) parameterises by `H_s`, `T_p` (ITTC/modified PM,
  JONSWAP, Torsethaugen). The code has no `H_s/T_p` interface; users must map
  sea state to `U10`. JONSWAP (peaked, fetch-limited) is absent.
* Directional spreading exists. The cos-2s normalisation
  `C_s = Gamma(s+1) / (2 sqrt(pi) Gamma(s+1/2))` equals Fossen's
  `2^{2s-1} Gamma^2(s+1) / (pi Gamma(2s+1))` by the duplication formula
  (verified). The ECKV spreading `(1 + Delta cos 2 phi)/(2 pi)` is front-back
  symmetric (standing-wave energy) and is only used when
  `use_symmetric_spreading_fn_` is set (default off). Correct.
* The trochoid algorithm's `angle` parameter spreads `number` components at
  `n * angle` from the mean direction: a crude discrete spreading, fine for
  its purpose.

---

## 6. Restoring `g(eta)` (Fossen Ch. 4) - verified (F-16)

Fossen Sec. 4.2 linearises to `G = diag(0, 0, rho g A_wp, rho g V GM_T, rho g V GM_L, 0)`.
The code integrates the hydrostatic pressure over the exact wetted triangles
(with the Kerner centre-of-pressure formulas `t_c = (4 z0 + 3h)/(6 z0 + 4h)` apex-up,
`(2 z0 + h)/(6 z0 + 2h)` apex-down). Using the port in `hydro_reference.py`:

| Check (box 4 x 1.5 x 0.6 m, 1000 kg, KG = 0.3 m) | code | `rho g V GM phi` |
|---|---|---|
| Heave: `F_z` at equilibrium draft | 9810.0 N | 9810.0 N |
| Roll restoring moment at `phi = 0.02 rad` (`GM_T = 0.934 m`) | -183.37 N m | -183.33 N m |
| Pitch restoring moment at `theta = 0.01 rad` (`GM_L = 8.0 m`) | -782.99 N m | -782.97 N m |
| Fully submerged rotated cube: force / moment about CoB | `rho g V` / `0` | `rho g V` / `0` |

The 0.02 % differences are the `sin phi ~ phi` linearisation, so the exact
nonlinear `g(eta)` is correct, `A_wp`, `GM_T` and `GM_L` emerge from the mesh,
and there is nothing to fix. Note that with `M_A = 0` (F-01) the *stiffness* is
right but the natural periods `2 pi sqrt((M + A)/C)` are not.

---

## 7. Numerics (F-11)

* Every hydrodynamic force is evaluated from the state at the start of the step
  and applied to DART as an external force; the code comments note that the
  link velocity components are "one time-step behind". The hydrodynamics are
  therefore an **explicit scheme with a one-step velocity delay** on top of
  DART's semi-implicit rigid-body integrator, with `max_step_size = 0.001 s` in
  the shipped worlds.
* Stability of the delayed explicit damping update
  `v_{n+1} = v_n - a v_{n-1}`, `a = dt d_eff / m`: roots of `z^2 - z + a` are inside
  the unit circle iff `0 < a < 1`. Hence the constraint per DOF
  ```
  dt * (c1_i + 2 c2_i |v_i| + d_tri_i(|v|)) / m_i < 1
  ```
  where `d_tri_i` is the linearised per-triangle drag slope (pressure drag
  `~ (cPDrag1 + 2 cPDrag2 v / vRDrag) S / vRDrag`, summed over faces). Without
  the delay the bound would be 2. For the WAM-V defaults the margin is large;
  it shrinks with the randomised `cSDrag2` up to 800 and with small light
  bodies (`m_i` small).
* Hydrostatic stiffness with a delayed force: `omega_n dt < ~1`
  (`omega_n = sqrt(rho g A_wp / m)`); at 1 ms this is not binding for any of
  the shipped models.
* Added mass: if `M_A` is introduced via `<fluid_added_mass>` (F-01 option 1)
  it becomes part of the engine's mass matrix and is integrated implicitly; no
  stiffness issue. If instead `-M_A nu_dot` were applied as an external force
  from a finite-difference `nu_dot`, the update `v_{n+1} = v_n + dt M_RB^-1 (-M_A (v_n - v_{n-1})/dt)`
  has amplification `M_RB^-1 M_A`, which is unstable for `M_A > M_RB` in any
  DOF. Do not do this.
* Foil lift with the current sign (F-02) is not dissipative-consistent with the
  pressure drag on the same face; on the faces it acts on it adds a normal
  force in the direction of the relative normal velocity, which reduces the
  effective `d_eff` above and can make the delayed explicit update diverge
  ("blow-over").

---

## 8. Minor findings

* **F-06 (moderate, applied)** Aero pass in `Hydrodynamics.cc`: selects faces
  with `vn = (v_wind - v_hull) . n > 0` and pushes them along `+n`. With outward
  normals this is the *lee* side. The net force on a symmetric body is the
  same, but the load is applied on the wrong faces, so the moment arm and the
  loaded area for asymmetric superstructures are wrong. Fossen Sec. 8.1 (wind
  forces) uses the projected *windward* area. Fixed to
  `vn = (v_hull - v_wind) . n > 0`, `F = -1/2 rho_air C_d A vn^2 n`.
* **F-12 (cosmetic)** `C_L_alpha = 2 pi cLift1` is the 2-D thin-foil slope; a
  low-aspect-ratio planing surface has `C_L_alpha = 2 pi / (1 + 2/AR)`
  (lifting-line, cf. Fossen's rudder lift in Ch. 9). `AR = B_wl^2 / A_bottom`
  uses the full waterline beam; for a catamaran that spans both hulls and
  over-predicts `AR` (under-predicts induced drag) by roughly `(B/b_hull)^2`.
* **F-13 (cosmetic)** `ComputeDampingForce` rotates into the body frame with
  `pose.Rot()` of the *CoM* pose (`linkPose * inertial.Pose()`), so a rotated
  inertial frame misaligns surge/sway/heave; use the link rotation.
* **F-14 (cosmetic)** `cLMax` default `cLift1 2 pi sin(alphaStall)` vs the linear
  branch `cLift1 2 pi alphaStall`: 1.1 % jump at stall.
* `ViscousDragCoefficient` (ITTC-57) and `DeepWaterDispersion*` (`omega^2 = g k`)
  are correct.
* `SetFromSDF` reads `foil_lift_on`, `cLift1`, `alphaStall`, `cLMax` *before*
  the `<randomize>` early-return, so they are honoured in randomised runs;
  fine.

---

## 9. What was applied (Phase 2)

| Commit | Finding | Behaviour change |
|---|---|---|
| `hydro: CoM velocity - Fossen Sec. 3.1 rigid-body kinematics` | F-05 | none when the inertial pose is at the link origin; otherwise drag/lift/damping torques change by the `omega x r_g` term |
| `hydro: wind load on windward faces - Fossen Sec. 8.1` | F-06 | only with `<aerodynamic_drag_on>`; net force unchanged for symmetric hulls, moments and asymmetric-hull magnitudes change |

Everything else is listed in the PR under "Not applied - needs decision".
