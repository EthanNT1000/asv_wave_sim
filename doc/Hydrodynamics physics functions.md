# Hydrodynamics Physics

[Hydrodynamics Physics](https://github.com/EthanNT1000/asv_wave_sim/blob/ament_environment_hooks/gz-waves/src/Physics.cc#L811)

1. Compute SubmergedTriangles, Areas, WaterlineLength
1. Compute point velocity at a triangles centroid

```c++
void Hydrodynamics::ComputePointVelocities()
{
  auto& position = this->data->position;
  auto& v = this->data->linVelocity;
  auto& omega = this->data->angVelocity;

  for (auto&& subTriProps : this->data->submergedTriangleProperties)
  {
    // relative position of the centroid wrt CoM
    subTriProps.xr = subTriProps.centroid - position;

    // vp = v + omega x xr
    subTriProps.vp = v + CGAL::cross_product(omega, subTriProps.xr);

    // up = vp / ||vp||
    subTriProps.up = Geometry::Normalize(subTriProps.vp);

    // cos(theta) = up . n
    subTriProps.cosTheta = CGAL::scalar_product(
        subTriProps.up, subTriProps.normal);

    // vn = (up . n) n
    subTriProps.vn = subTriProps.normal * subTriProps.cosTheta;

    // vt = vp - vn
    subTriProps.vt = subTriProps.vp - subTriProps.vn;

    // un = vn / ||vn||
    // subTriProps.un = Geometry::Normalize(subTriProps.vn);

    // ut = vt / ||vt||
    subTriProps.ut = Geometry::Normalize(subTriProps.vt);

    // uf = - vt / ||vt|| = - ut
    subTriProps.uf = - subTriProps.ut;

    // vf = ||vp|| uf
    subTriProps.vf =
        subTriProps.uf * std::sqrt(subTriProps.vp.squared_length());
  }
}
```

1. Compute Buoyancy Force

```c++
void Hydrodynamics::ComputeBuoyancyForce()
{
  cgal::Vector3 sumForce  = CGAL::NULL_VECTOR;
  cgal::Vector3 sumTorque = CGAL::NULL_VECTOR;

  this->data->fBuoyancy.clear();
  this->data->cBuoyancy.clear();

  // Calculate the buoyancy force for the submerged triangles
  auto& position  =  this->data->position;
  auto& wavefieldSampler = *this->data->wavefieldSampler;
  for (auto&& subTri : this->data->submergedTriangles)
  {
    // Force and center of pressure.
    cgal::Point3 center = CGAL::ORIGIN;
    cgal::Vector3 force = CGAL::NULL_VECTOR;
    Physics::BuoyancyForceAtCenterOfPressure(
      wavefieldSampler, subTri, center, force);
    this->data->fBuoyancy.push_back(force);
    this->data->cBuoyancy.push_back(center);

    // Torque
    cgal::Vector3 xr = center - position;
    cgal::Vector3 torque = CGAL::cross_product(xr, force);
    sumForce += force;
    sumTorque += torque;
  }

  this->data->force  += sumForce;
  this->data->torque += sumTorque;
}

void Physics::BuoyancyForceAtCenterOfPressure(
  const WavefieldSampler& _wavefieldSampler,
  const cgal::Triangle& _triangle,
  cgal::Point3& _center,
  cgal::Vector3& _force
)
{
  // Sort triangle vertices by height.
  std::array<cgal::Point3, 3> v {
    _triangle[0],
    _triangle[1],
    _triangle[2]
  };
  std::array<double, 3> vz { v[0].z(), v[1].z(), v[2].z() };
  auto index = algorithm::sort_indexes(vz);

  cgal::Point3 H = v[index[0]];
  cgal::Point3 M = v[index[1]];
  cgal::Point3 L = v[index[2]];

  // @DEBUG_INFO
  // DebugPrint(_triangle);
  // gzmsg << "vz:          "; for (auto z: vz)    { gzmsg << z << " "; };
  // gzmsg << "\n";
  // gzmsg << "index:       "; for (auto i: index) { gzmsg << i << " "; };
  // gzmsg << "\n";
  // gzmsg << "H:           " << H << "\n";
  // gzmsg << "M:           " << M << "\n";
  // gzmsg << "L:           " << L << "\n";

  // Calculate the depth at the centroid
  cgal::Point3 C = Geometry::TriangleCentroid(_triangle);

  // Calculate the depth
  double depthC = _wavefieldSampler.ComputeDepth(C);

  cgal::Vector3 normal = Geometry::Normal(_triangle);

  // Calculate buoyancy
  BuoyancyForceAtCenterOfPressure(depthC, C, H, M, L, normal, _center, _force);
}

void Physics::BuoyancyForceAtCenterOfPressure(
  double _depthC,
  const cgal::Point3& _C,
  const cgal::Point3& _H,
  const cgal::Point3& _M,
  const cgal::Point3& _L,
  const cgal::Vector3& _normal,
  cgal::Point3& _center,
  cgal::Vector3& _force
)
{
  double fluidDensity = PhysicalConstants::WaterDensity();  // kg m^-3
  double gravity = PhysicalConstants::Gravity();            // m s^-1

  // Split the triangle into upper and lower triangles bisected by a
  // line normal to the z-axis
  cgal::Point3 D = Geometry::HorizontalIntercept(_H, _M, _L);
  cgal::Point3 B = Geometry::MidPoint(_M, D);

  // Initialise to the base midpoint (correct force calcuation for
  // triangles with a horizontal base)
  double fU = 0, fL = 0;
  cgal::Point3 CpU = B;
  cgal::Point3 CpL = B;

  // Upper triangle H > M
  if (_H.z() >= _M.z())
  {
    // Center of pressure
    double z0 = _depthC - (_H.z() - _C.z());
    CpU = CenterOfPressureApexUp(z0, _H, _M, B);

    // Force at centroid
    cgal::Point3 CU = Geometry::TriangleCentroid(_H, _M, D);
    double hCU = _depthC + (_C.z() - CU.z());
    double area = Geometry::TriangleArea(_H, _M, D);
    fU = fluidDensity * gravity * area * hCU;

    // @DEBUG_INFO
    // gzmsg << "_depthC: " << _depthC << "\n";
    // gzmsg << "_C:      " << _C << "\n";
    // gzmsg << "_H:      " << _H << "\n";
    // gzmsg << "_M:      " << _M << "\n";
    // gzmsg << "_L:      " << _L << "\n";
    // gzmsg << "_normal: " << _normal << "\n";
    // gzmsg << "D:       " << D << "\n";
    // gzmsg << "B:       " << B << "\n";
    // gzmsg << "z0:      " << z0 << "\n";
    // gzmsg << "CpU:     " << CpU << "\n";
    // gzmsg << "CU:      " << CU << "\n";
    // gzmsg << "hCU:     " << hCU << "\n";
    // gzmsg << "area:    " << area << "\n";
    // gzmsg << "fU:      " << fU << "\n";
  }

  // Lower triangle L < M
  if (_M.z() > _L.z())
  {
    // Center of preseesure
    double z0 = _depthC + (_C.z() - _M.z());
    CpL = CenterOfPressureApexDn(z0, _L, _M, B);

    // Force at centroid
    cgal::Point3 CL = Geometry::TriangleCentroid(_L, _M, D);
    double hCL = _depthC + (_C.z() - CL.z());
    double area = Geometry::TriangleArea(_L, _M, D);
    fL = fluidDensity * gravity * area * hCL;

    // @DEBUG_INFO
    // gzmsg << "_depthC: " << _depthC << "\n";
    // gzmsg << "_C:      " << _C << "\n";
    // gzmsg << "_H:      " << _H << "\n";
    // gzmsg << "_M:      " << _M << "\n";
    // gzmsg << "_L:      " << _L << "\n";
    // gzmsg << "_normal: " << _normal << "\n";
    // gzmsg << "D:       " << D << "\n";
    // gzmsg << "B:       " << B << "\n";
    // gzmsg << "z0:      " << z0 << "\n";
    // gzmsg << "CpL:     " << CpL << "\n";
    // gzmsg << "CL:      " << CL << "\n";
    // gzmsg << "hCL:     " << hCL << "\n";
    // gzmsg << "area:    " << area << "\n";
    // gzmsg << "fL:      " << fL << "\n";
  }

  // Calculate the force and centre of application
  _force = _normal * (fU + fL);
  _center = CenterOfForce(fU, fL, CpU, CpL);

  // @DEBUG_INFO
  // gzmsg << "_depthC: " << _depthC << "\n";
  // gzmsg << "_C:      " << _C << "\n";
  // gzmsg << "_H:      " << _H << "\n";
  // gzmsg << "_M:      " << _M << "\n";
  // gzmsg << "_L:      " << _L << "\n";
  // gzmsg << "_normal: " << _normal << "\n";
  // gzmsg << "_force:  " << _force << "\n";
  // gzmsg << "_center: " << _center << "\n";
}

cgal::Point3 Physics::CenterOfPressureApexUp(
  double _z0,
  const cgal::Point3& _H,
  const cgal::Point3& _M,
  const cgal::Point3& _B
)
{
  cgal::Vector3 alt = _B - _H;
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

//////////////////////////////////////////////////
cgal::Point3 Physics::CenterOfPressureApexDn(
  double _z0,
  const cgal::Point3& _L,
  const cgal::Point3& _M,
  const cgal::Point3& _B
)
{
  cgal::Vector3 alt = _L - _B;
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

cgal::Point3 Physics::CenterOfForce(
  double _fA, double _fB,
  const cgal::Point3& _A,
  const cgal::Point3& _B
)
{
  /// \todo provide robust floating point checks
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

1. Compute Viscous Drag Force

```c++
void Hydrodynamics::ComputeViscousDragForce()
{
  double rho = PhysicalConstants::WaterDensity();
  double Rn = this->ComputeReynoldsNumber();
  double cF = Physics::ViscousDragCoefficient(Rn);

  cgal::Vector3 sumForce  = CGAL::NULL_VECTOR;
  cgal::Vector3 sumTorque = CGAL::NULL_VECTOR;
  for (auto&& subTriProps : this->data->submergedTriangleProperties)
  {
    // Force
    double fDrag = 0.5 * rho * cF * subTriProps.area
      * std::sqrt(subTriProps.vf.squared_length());
    cgal::Vector3 force = subTriProps.vf * fDrag;
    sumForce += force;

    // Torque;
    cgal::Vector3 torque = CGAL::cross_product(subTriProps.xr, force);
    sumTorque += torque;
  }

  this->data->force  += sumForce;
  this->data->torque += sumTorque;
}

double Hydrodynamics::ComputeReynoldsNumber() const
{
  // fluid speed
  auto& v = this->data->linVelocity;
  double u = std::sqrt(v.squared_length());

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
```

1. Compute Pressure Drag Force

```c++
void Hydrodynamics::ComputePressureDragForce()
{
  auto& params = *this->data->params;

  // Positive pressure
  double cPDrag1 = params.CPDrag1();
  double cPDrag2 = params.CPDrag2();
  double fPDrag  = params.FPDrag();
  // Negative pressure (suction)
  double cSDrag1 = params.CSDrag1();
  double cSDrag2 = params.CSDrag2();
  double fSDrag  = params.FSDrag();

  // Reference speed
  double vRDrag  = params.VRDrag();

  cgal::Vector3 sumForce  = CGAL::NULL_VECTOR;
  cgal::Vector3 sumTorque = CGAL::NULL_VECTOR;
  for (auto&& subTriProps : this->data->submergedTriangleProperties)
  {
    // General
    double S    = subTriProps.area;
    double vp   = std::sqrt(subTriProps.vp.squared_length());
    double cosTheta = subTriProps.cosTheta;

    double v    = vp / vRDrag;
    double drag = 0.0;
    if (cosTheta >= 0.0)
    {
      drag = -(cPDrag1 * v + cPDrag2 * v * v) * S * std::pow(cosTheta, fPDrag);
    } else {
      drag =  (cSDrag1 * v + cSDrag2 * v * v) * S * std::pow(-cosTheta, fSDrag);
    }
    cgal::Vector3 force = subTriProps.normal * drag;
    sumForce += force;

    // Torque;
    cgal::Vector3 torque = CGAL::cross_product(subTriProps.xr, force);
    sumTorque += torque;
  }

  this->data->force  += sumForce;
  this->data->torque += sumTorque;

  // @DEBUG_INFO
  // if (std::abs(sumForce.z()) > 1.0E+10)
  // {
  //   gzmsg << "Overflow in ComputePressureDragForce..."    << "\n";
  //   gzmsg << "position:     " << this->data->position     << "\n";
  //   gzmsg << "linVelocity:  " << this->data->linVelocity  << "\n";
  //   gzmsg << "angVelocity:  " << this->data->angVelocity  << "\n";
  //   gzmsg << "force:        " << sumForce                 << "\n";
  //   gzmsg << "torque:       " << sumTorque                << "\n";
  //   for (auto&& subTriProps : this->data->submergedTriangleProperties)
  //   {
  //     DebugPrint(subTriProps);
  //   }
  // }
}
```

1. Compute Damping Force

```c++
void Hydrodynamics::ComputeDampingForce()
{
  auto& params = *this->data->params;

    // Linear drag coefficients
  double cDampL1 = params.CDampL1();
  double cDampR1 = params.CDampR1();

  // Quadratic drag coefficients
  double cDampL2 = params.CDampL2();
  double cDampR2 = params.CDampR2();

  double area = this->data->area;
  double subArea = this->data->submergedArea;
  double rs = subArea / area;

  // Force
  auto& v = this->data->linVelocity;
  double linSpeed = std::sqrt(v.squared_length());
  double cL = - rs * (cDampL1 + cDampL2 * linSpeed);
  cgal::Vector3 force = v * cL;

  auto& omega = this->data->angVelocity;
  double angSpeed = std::sqrt(omega.squared_length());
  double cR = - rs * (cDampR1 + cDampR2 * angSpeed);
  cgal::Vector3 torque = omega * cR;

  this->data->force  += force;
  this->data->torque += torque;

  // @DEBUG_INF0
  // gzmsg << "area:       " << area << "\n";
  // gzmsg << "subArea:    " << subArea << "\n";
  // gzmsg << "v:          " << v << "\n";
  // gzmsg << "omega:      " << omega << "\n";
  // gzmsg << "linSpeed:   " << linSpeed << "\n";
  // gzmsg << "angSpeed:   " << angSpeed << "\n";
  // gzmsg << "cR:         " << cR << "\n";
  // gzmsg << "cL:         " << cL << "\n";
  // gzmsg << "force:      " << force << "\n";
  // gzmsg << "torque:     " << torque << "\n";
}
```