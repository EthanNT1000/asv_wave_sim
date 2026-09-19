# Hydrodynamics Physics

[Hydrodynamics Physics](https://github.com/EthanNT1000/asv_wave_sim/blob/ament_environment_hooks/gz-waves/src/Physics.cc#L838)

See [`Hydrodynamics Physics Algorithm.md`](Hydrodynamics%20Physics%20Algorithm.md) for the derivation behind each of these functions. This file mirrors the current source in `gz-waves/src/Physics.cc` (per-triangle force loop) and `gz-waves/src/systems/hydrodynamics/Hydrodynamics.cc` (above-waterline aerodynamic drag).

1. Compute SubmergedTriangles, Areas, WaterlineLength, WaterlineBeam, dynamic foil geometry, and sample the bulk current at the CoM — see `Hydrodynamics::Update`.
1. Single parallel pass over submerged triangles computing point velocities, buoyancy, viscous drag, pressure drag and foil lift:

```c++
void Hydrodynamics::ComputeAllSubmergedForces(
    const std::chrono::_V2::steady_clock::duration& simTime)
{
  const int n = static_cast<int>(this->data->submergedTriangleProperties.size());
  if (n == 0) return;

  this->data->fBuoyancy.resize(n);
  this->data->cBuoyancy.resize(n);

  auto& position         = this->data->position;
  auto& v_body           = this->data->linVelocity;
  auto& omega            = this->data->angVelocity;
  auto& wavefieldSampler = *this->data->wavefieldSampler;
  const double t         = std::chrono::duration<double>(simTime).count();

  // Viscous drag
  const double rho  = PhysicalConstants::WaterDensity();
  const double Rn   = this->ComputeReynoldsNumber();
  const double cF   = Physics::ViscousDragCoefficient(Rn);
  const bool viscousOn  = this->data->params->ViscousDragOn();

  // Pressure drag
  const bool pressureOn = this->data->params->PressureDragOn();
  const double cPDrag1  = this->data->params->CPDrag1();
  const double cPDrag2  = this->data->params->CPDrag2();
  const double fPDrag   = this->data->params->FPDrag();
  const double cSDrag1  = this->data->params->CSDrag1();
  const double cSDrag2  = this->data->params->CSDrag2();
  const double fSDrag   = this->data->params->FSDrag();
  const double vRDrag   = this->data->params->VRDrag();

  // Foil lift
  const bool foilOn        = this->data->params->FoilLiftOn();
  const double Cl_alpha    = this->data->params->CLift1() * 2.0 * M_PI;
  const double alpha_stall = this->data->params->AlphaStall();
  const double Cl_max      = this->data->params->CLMax();
  double AR = std::max(this->data->dynamic_foil_ar, 0.5);
  const double bottomThresh = this->data->params->BOTTOM_THRESHOLD;

  double fx = 0.0, fy = 0.0, fz = 0.0;
  double tx = 0.0, ty = 0.0, tz = 0.0;

#pragma omp parallel for reduction(+:fx,fy,fz,tx,ty,tz) schedule(static)
  for (int i = 0; i < n; ++i)
  {
    auto& props          = this->data->submergedTriangleProperties[i];
    const auto& subTri   = this->data->submergedTriangles[i];

    // ── Point velocities ─────────────────────────────────────────────────
    ComputePointVelocities(props, position, v_body, omega,
        wavefieldSampler, t, this->data->params->GetWaterCurrentGrid());

    // ── Helper: accumulate a (force, torque) pair into the reduction scalars ──
    auto acc = [&](const std::pair<geom::Vector3, geom::Vector3>& ft) {
      fx += geom::ToDouble(ft.first.x());
      fy += geom::ToDouble(ft.first.y());
      fz += geom::ToDouble(ft.first.z());
      tx += geom::ToDouble(ft.second.x());
      ty += geom::ToDouble(ft.second.y());
      tz += geom::ToDouble(ft.second.z());
    };

    // ── Buoyancy ──────────────────────────────────────────────────────────
    acc(ComputeBuoyancyForce(
        wavefieldSampler, subTri, position,
        this->data->fBuoyancy[i], this->data->cBuoyancy[i]));

    // ── Viscous drag ──────────────────────────────────────────────────────
    if (viscousOn)
      acc(ComputeViscousDragForce(props, rho, cF));

    // ── Pressure drag ─────────────────────────────────────────────────────
    if (pressureOn)
      acc(ComputePressureDragForce(
          props, cPDrag1, cPDrag2, fPDrag, cSDrag1, cSDrag2, fSDrag, vRDrag));

    // ── Foil lift ─────────────────────────────────────────────────────────
    if (foilOn)
      acc(ComputeFoilLiftForce(
          props, rho, Cl_alpha, alpha_stall, Cl_max, AR, bottomThresh));
  }

  this->data->force  += geom::Vector3(fx, fy, fz);
  this->data->torque += geom::Vector3(tx, ty, tz);
}
```

1. Compute point velocity and the hull-relative-to-fluid velocity at a triangle's centroid — note this is now a `static` function taking all its inputs as parameters (called from inside the parallel loop above), rather than a member function reading directly from `this->data`:

```c++
void Hydrodynamics::ComputePointVelocities(
    SubmergedTriangleProperties& props,
    const geom::Point3& position,
    const geom::Vector3& v_body,
    const geom::Vector3& omega,
    const WavefieldSampler& wavefieldSampler,
    double t,
    const WaterCurrentGrid& currentGrid)
{
  props.xr = props.centroid - position;
  props.vp = v_body + geom::Cross(omega, props.xr);

  const double cx = props.centroid.x();
  const double cy = props.centroid.y();
  const double cz = props.centroid.z();

  props.v_orbital = wavefieldSampler.ComputeOrbitalVelocity(cx, cy, cz, t);
  props.v_current = currentGrid.SampleAt(cx, cy);
  props.v_fluid   = props.v_orbital + props.v_current;
  props.v_rel     = props.vp - props.v_fluid;
  props.v_rel_mag = std::sqrt(geom::ToDouble(props.v_rel.squared_length()));

  const double v_rel_dot_n = geom::ToDouble(
      geom::Dot(props.v_rel, props.normal));
  props.v_rel_n = props.normal * v_rel_dot_n;
  props.v_rel_t = props.v_rel - props.v_rel_n;

  const double v_rel_t_mag = std::sqrt(
      geom::ToDouble(props.v_rel_t.squared_length()));
  props.alpha = std::atan2(v_rel_dot_n, v_rel_t_mag + 1e-9);

  props.up = (props.v_rel_mag > 1e-9)
      ? props.v_rel / props.v_rel_mag
      : geom::NullVector();
  props.cosTheta = geom::Dot(props.up, props.normal);
  props.vn = props.normal * props.cosTheta * props.v_rel_mag;
  props.vt = props.v_rel - props.vn;
  props.ut = (v_rel_t_mag > 1e-9)
      ? props.v_rel_t / v_rel_t_mag
      : geom::NullVector();
  props.uf = -props.ut;
  props.vf = props.uf * props.v_rel_mag;
}
```

1. Compute Buoyancy Force — the outer wrapper now returns a `(force, torque)` pair for the caller to accumulate; the inner center-of-pressure math (`Physics::BuoyancyForceAtCenterOfPressure` and its helpers) is unchanged from the original model:

```c++
std::pair<geom::Vector3, geom::Vector3>
Hydrodynamics::ComputeBuoyancyForce(
    const WavefieldSampler& wavefieldSampler,
    const geom::Triangle& subTri,
    const geom::Point3& position,
    geom::Vector3& bForce_out,
    geom::Point3& bCenter_out)
{
  geom::Point3  bCenter = geom::Origin();
  geom::Vector3 bForce  = geom::NullVector();
  Physics::BuoyancyForceAtCenterOfPressure(
      wavefieldSampler, subTri, bCenter, bForce);
  bForce_out  = bForce;
  bCenter_out = bCenter;
  const geom::Vector3 xr     = bCenter - position;
  const geom::Vector3 torque = geom::Cross(xr, bForce);
  return {bForce, torque};
}

void Physics::BuoyancyForceAtCenterOfPressure(
  const WavefieldSampler& _wavefieldSampler,
  const geom::Triangle& _triangle,
  geom::Point3& _center,
  geom::Vector3& _force
)
{
  // Sort triangle vertices by height.
  std::array<geom::Point3, 3> v {
    _triangle[0],
    _triangle[1],
    _triangle[2]
  };
  std::array<double, 3> vz { v[0].z(), v[1].z(), v[2].z() };
  auto index = algorithm::sort_indexes(vz);

  geom::Point3 H = v[index[0]];
  geom::Point3 M = v[index[1]];
  geom::Point3 L = v[index[2]];

  // Calculate the depth at the centroid
  geom::Point3 C = Geometry::TriangleCentroid(_triangle);
  double depthC = _wavefieldSampler.ComputeDepth(C);
  geom::Vector3 normal = Geometry::Normal(_triangle);

  // Calculate buoyancy
  BuoyancyForceAtCenterOfPressure(depthC, C, H, M, L, normal, _center, _force);
}

void Physics::BuoyancyForceAtCenterOfPressure(
  double _depthC,
  const geom::Point3& _C,
  const geom::Point3& _H,
  const geom::Point3& _M,
  const geom::Point3& _L,
  const geom::Vector3& _normal,
  geom::Point3& _center,
  geom::Vector3& _force
)
{
  double fluidDensity = PhysicalConstants::WaterDensity();  // kg m^-3
  double gravity = PhysicalConstants::Gravity();            // m s^-1

  // Split the triangle into upper and lower triangles bisected by a
  // line normal to the z-axis
  geom::Point3 D = Geometry::HorizontalIntercept(_H, _M, _L);
  geom::Point3 B = Geometry::MidPoint(_M, D);

  double fU = 0, fL = 0;
  geom::Point3 CpU = B;
  geom::Point3 CpL = B;

  // Upper triangle H > M
  if (_H.z() >= _M.z())
  {
    double z0 = _depthC - (_H.z() - _C.z());
    CpU = CenterOfPressureApexUp(z0, _H, _M, B);

    geom::Point3 CU = Geometry::TriangleCentroid(_H, _M, D);
    double hCU = _depthC + (_C.z() - CU.z());
    double area = Geometry::TriangleArea(_H, _M, D);
    fU = fluidDensity * gravity * area * hCU;
  }

  // Lower triangle L < M
  if (_M.z() > _L.z())
  {
    double z0 = _depthC + (_C.z() - _M.z());
    CpL = CenterOfPressureApexDn(z0, _L, _M, B);

    geom::Point3 CL = Geometry::TriangleCentroid(_L, _M, D);
    double hCL = _depthC + (_C.z() - CL.z());
    double area = Geometry::TriangleArea(_L, _M, D);
    fL = fluidDensity * gravity * area * hCL;
  }

  // Calculate the force and centre of application
  _force = _normal * (fU + fL);
  _center = CenterOfForce(fU, fL, CpU, CpL);
}

geom::Point3 Physics::CenterOfPressureApexUp(
  double _z0,
  const geom::Point3& _H,
  const geom::Point3& _M,
  const geom::Point3& _B
)
{
  geom::Vector3 alt = _B - _H;
  double h = _H.z() - _M.z();
  double tc = 2.0/3.0;
  double div = 6.0 * _z0 + 4.0 * h;
  constexpr double tol = 1.0E-16;
  if (std::fabs(div) > tol)
  {
    tc = (4.0 * _z0 + 3.0 * h) / div;
  }
  return _H + alt * tc;
}

geom::Point3 Physics::CenterOfPressureApexDn(
  double _z0,
  const geom::Point3& _L,
  const geom::Point3& _M,
  const geom::Point3& _B
)
{
  geom::Vector3 alt = _L - _B;
  double h = _M.z() - _L.z();
  double tc = 1.0/3.0;
  double div = 6.0 * _z0 + 2.0 * h;
  constexpr double tol = 1.0E-16;
  if (std::fabs(div) > tol)
  {
    tc = (2.0 * _z0 + h) / div;
  }
  return _B + alt * tc;
}

geom::Point3 Physics::CenterOfForce(
  double _fA, double _fB,
  const geom::Point3& _A,
  const geom::Point3& _B
)
{
  double div = _fA + _fB;
  constexpr double tol = 1.0E-16;
  if (std::fabs(div) > tol)
  {
    double t = _fA / div;
    return _B + (_A - _B) * t;
  } else {
    return _B;
  }
}
```

1. Compute Viscous Drag Force — now driven by `props.v_rel_mag` (hull speed relative to wave orbital velocity + bulk current), returned as a `(force, torque)` pair instead of accumulated into member state:

```c++
double Hydrodynamics::ComputeReynoldsNumber() const
{
  // fluid speed, relative to the bulk current sampled at the CoM
  geom::Vector3 v_rel = this->data->linVelocity - this->data->waterCurrentCoM;
  double u = std::sqrt(v_rel.squared_length());

  // characteristic length
  double L = this->data->waterlineLength;

  // Reynolds number
  double kv = PhysicalConstants::WaterKinematicViscosity();
  double Rn = u * L / kv;

  return Rn;
}

double Physics::ViscousDragCoefficient(double Rn)
{
  // Set a lower limit on Rn to 1.0E+3 since the 1957 ITTC formula
  // has a pole at Rn = 1.0E+2. For a 10m boat in salt water this
  // corresponds to a velocity ~ 1.0E-5
  double r = std::max(1.0E3, Rn);
  double d = std::log10(r) - 2.0;
  double d2 = d*d;
  double CF = 0.075 / d2;
  return CF;
}

// Viscous drag force - applied at triangle centroid.
std::pair<geom::Vector3, geom::Vector3>
Hydrodynamics::ComputeViscousDragForce(
    const SubmergedTriangleProperties& props,
    double rho,
    double cF)
{
  const double fDrag         = 0.5 * rho * cF * props.area * props.v_rel_mag;
  const geom::Vector3 force  = props.vf * fDrag;
  const geom::Vector3 torque = geom::Cross(props.xr, force);
  return {force, torque};
}
```

1. Compute Pressure Drag Force — same functional form as before, but scaled by `props.v_rel_mag` instead of the raw point-velocity magnitude:

```c++
// Pressure drag force - applied at triangle centroid.
std::pair<geom::Vector3, geom::Vector3>
Hydrodynamics::ComputePressureDragForce(
    const SubmergedTriangleProperties& props,
    double cPDrag1, double cPDrag2, double fPDrag,
    double cSDrag1, double cSDrag2, double fSDrag,
    double vRDrag)
{
  const double S        = props.area;
  const double v        = props.v_rel_mag / vRDrag;
  const double cosTheta = geom::ToDouble(props.cosTheta);
  const double drag = (cosTheta >= 0.0)
      ? -(cPDrag1 * v + cPDrag2 * v * v) * S * std::pow(cosTheta,  fPDrag)
      :  (cSDrag1 * v + cSDrag2 * v * v) * S * std::pow(-cosTheta, fSDrag);
  const geom::Vector3 force  = props.normal * drag;
  const geom::Vector3 torque = geom::Cross(props.xr, force);
  return {force, torque};
}
```

1. Compute Foil Lift Force — new: dynamic lift on near-horizontal, bottom-facing submerged triangles, with a stall cap on $C_l$ and an induced-drag penalty based on the dynamic foil aspect ratio:

```c++
std::pair<geom::Vector3, geom::Vector3>
Hydrodynamics::ComputeFoilLiftForce(
    const SubmergedTriangleProperties& props,
    double rho,
    double Cl_alpha, double alpha_stall, double Cl_max,
    double AR, double bottomThresh)
{
  const double nz = geom::ToDouble(props.normal.z());
  if (nz < bottomThresh || props.v_rel_mag < 1e-4)
    return {geom::NullVector(), geom::NullVector()};

  const double alpha = props.alpha;
  const double Cl = (std::fabs(alpha) < alpha_stall)
      ? Cl_alpha * alpha
      : Cl_max * (alpha > 0.0 ? 1.0 : -1.0);
  const double Cdi = (Cl * Cl) / (M_PI * AR + 1e-9);
  const double q_A = 0.5 * rho * props.v_rel_mag * props.v_rel_mag * props.area;

  geom::Vector3 lift_dir = props.normal
      - props.up * geom::ToDouble(
            geom::Dot(props.normal, props.up));
  const double ld_mag = std::sqrt(
      geom::ToDouble(lift_dir.squared_length()));
  if (ld_mag < 1e-9)
    return {geom::NullVector(), geom::NullVector()};

  lift_dir = lift_dir / ld_mag;
  const geom::Vector3 F_foil =
      lift_dir * (Cl * q_A) + (-props.up) * (Cdi * q_A);
  return {F_foil, geom::Cross(props.xr, F_foil)};
}
```

1. Compute Damping Force — replaced the single scalar linear/angular model with a full per-DOF (Fossen-style) damping matrix applied in the body frame, computed once per link per step (not per submerged triangle, and no longer scaled by the submerged-area ratio):

```c++
void Hydrodynamics::ComputeDampingForce()
{
  auto& params = *this->data->params;

  // Transform velocities to body frame for per-DOF Fossen damping matrix.
  // The gz::math quaternion is the world-to-body rotation stored in pose.Rot().
  gz::math::Quaterniond R     = this->data->pose.Rot();
  gz::math::Quaterniond R_inv = R.Inverse();

  gz::math::Vector3d linW = ToGz(this->data->linVelocity
                                  - this->data->waterCurrentCoM);
  gz::math::Vector3d angW = ToGz(this->data->angVelocity);

  gz::math::Vector3d linB = R_inv.RotateVector(linW);
  gz::math::Vector3d angB = R_inv.RotateVector(angW);

  // F = -(D1·v + D2·|v|·v) applied independently per DOF.
  auto damp = [](double c1, double c2, double v) -> double {
    return -(c1 * v + c2 * std::fabs(v) * v);
  };

  gz::math::Vector3d forceB(
    damp(params.CDampU1(), params.CDampU2(), linB.X()),  // surge
    damp(params.CDampV1(), params.CDampV2(), linB.Y()),  // sway
    damp(params.CDampW1(), params.CDampW2(), linB.Z())   // heave
  );
  gz::math::Vector3d torqueB(
    damp(params.CDampP1(), params.CDampP2(), angB.X()),  // roll
    damp(params.CDampQ1(), params.CDampQ2(), angB.Y()),  // pitch
    damp(params.CDampN1(), params.CDampN2(), angB.Z())   // yaw
  );

  this->data->force  += ToVector3(R.RotateVector(forceB));
  this->data->torque += ToVector3(R.RotateVector(torqueB));
}
```

1. Compute above-waterline aerodynamic drag — new, lives in the Gazebo system (`gz-waves/src/systems/hydrodynamics/Hydrodynamics.cc`), not `Physics.cc`. Runs over each hull's already-computed `TriangleProperties`, using the world `Wind` entity's velocity, in its own OpenMP reduction loop:

```c++
if (this->aeroDragOn)
{
  static constexpr double kRhoAir = 1.225;  // kg/m^3

  geom::Vector3 coMVec = waves::ToVector3(linkCoMPose.Pos());

  const auto& tris = hd->hydrodynamics[j]->GetTriangleProperties();
  const int nTris = static_cast<int>(tris.size());

  double fx = 0.0, fy = 0.0, fz = 0.0;
  double tx = 0.0, ty = 0.0, tz = 0.0;

  #pragma omp parallel for reduction(+:fx,fy,fz,tx,ty,tz) schedule(static)
  for (int ti = 0; ti < nTris; ++ti)
  {
    const auto& tri = tris[ti];

    // Determine above-waterline area. hh, hm, hl are vertex heights above
    // the water surface (sorted). Fully above-water faces leave subArea
    // as NaN (PopulateSubmergedTriangle returns early), so each case is
    // handled explicitly.
    double areaAbove;
    if      (tri.hh <= 0.0) continue;                          // fully submerged
    else if (tri.hl  > 0.0) areaAbove = tri.area;              // fully above
    else                    areaAbove = tri.area - tri.subArea; // partial

    if (areaAbove < 1e-9) continue;

    double nLen = std::sqrt(geom::ToDouble(tri.normal.squared_length()));
    if (nLen < 1e-9) continue;
    geom::Vector3 nHat = tri.normal / nLen;

    // Hull velocity at the triangle centroid (CoM + angular correction).
    geom::Vector3 centroid = geom::Vector3(
        (tri.vh.x() + tri.vm.x() + tri.vl.x()) / 3.0,
        (tri.vh.y() + tri.vm.y() + tri.vl.y()) / 3.0,
        (tri.vh.z() + tri.vm.z() + tri.vl.z()) / 3.0);
    geom::Vector3 r = centroid - coMVec;
    geom::Vector3 vHull = linVelocity + geom::Cross(angVelocity, r);

    // Normal component of relative wind.
    geom::Vector3 vRel = windVelocity - vHull;
    double vn = geom::ToDouble(vRel * nHat);
    if (vn <= 0.0) continue;  // lee side — no pressure

    geom::Vector3 f   = (0.5 * kRhoAir * this->cAeroDrag * areaAbove * vn * vn) * nHat;
    geom::Vector3 tau = geom::Cross(r, f);

    fx += geom::ToDouble(f.x()); fy += geom::ToDouble(f.y()); fz += geom::ToDouble(f.z());
    tx += geom::ToDouble(tau.x()); ty += geom::ToDouble(tau.y()); tz += geom::ToDouble(tau.z());
  }

  gz::math::Vector3d aeroForceSum(fx, fy, fz);
  gz::math::Vector3d aeroTorqueSum(tx, ty, tz);
  if (aeroForceSum.IsFinite())  hd->link.AddWorldForce(_ecm, aeroForceSum);
  if (aeroTorqueSum.IsFinite()) hd->link.AddWorldWrench(_ecm, gz::math::Vector3d::Zero, aeroTorqueSum);
}
```

1. Bilinear sampling of the bulk water-current grid (`gz-waves/src/WaterCurrentGrid.cc`), used by both `ComputePointVelocities` (per submerged triangle) and `SampleWaterCurrentCoM` (once per link, for Reynolds number and per-DOF damping):

```c++
gz::geom::Vector3 WaterCurrentGrid::SampleAt(double x, double y) const
{
  if (!loaded_) return {0.0, 0.0, 0.0};

  // Map world coords to fractional grid indices
  double fx = (x - x_min_) / cell_size_x_;
  double fy = (y - y_min_) / cell_size_y_;

  // Out-of-bounds -> zero current (no extrapolation)
  if (fx < 0.0 || fy < 0.0 || fx > nx_ - 1 || fy > ny_ - 1)
    return {0.0, 0.0, 0.0};

  // Bilinear interpolation
  int ix0 = static_cast<int>(fx);
  int iy0 = static_cast<int>(fy);
  int ix1 = std::min(ix0 + 1, nx_ - 1);
  int iy1 = std::min(iy0 + 1, ny_ - 1);

  double tx = fx - ix0;   // [0, 1)
  double ty = fy - iy0;

  auto cell = [&](int ix, int iy) -> std::pair<double, double> {
    const size_t idx = (static_cast<size_t>(iy) * nx_ + ix) * 2;
    return { data_[idx], data_[idx + 1] };
  };

  auto [vx00, vy00] = cell(ix0, iy0);
  auto [vx10, vy10] = cell(ix1, iy0);
  auto [vx01, vy01] = cell(ix0, iy1);
  auto [vx11, vy11] = cell(ix1, iy1);

  double vx = (1-tx)*(1-ty)*vx00 + tx*(1-ty)*vx10 + (1-tx)*ty*vx01 + tx*ty*vx11;
  double vy = (1-tx)*(1-ty)*vy00 + tx*(1-ty)*vy10 + (1-tx)*ty*vy01 + tx*ty*vy11;

  return { vx, vy, 0.0 };
}
```
