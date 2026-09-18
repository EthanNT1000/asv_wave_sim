// Copyright (C) 2019  Rhys Mainwaring
//
// This program is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with this program.  If not, see <https://www.gnu.org/licenses/>.

#include "gz/waves/Physics.hh"

#include <algorithm>
#include <array>
#include <cmath>
#include <iostream>
#include <iterator>
#include <limits>
#include <memory>
#include <string>
#include <utility>
#include <vector>

#include <gz/common/Console.hh>
#include <gz/math/Pose3.hh>
#include <gz/math/Vector3.hh>

#include <sdf/sdf.hh>

#include "gz/waves/Algorithm.hh"
#include "gz/waves/Convert.hh"
#include "gz/waves/Geometry.hh"
#include "gz/waves/PhysicalConstants.hh"
#include "gz/waves/Utilities.hh"
#include "gz/waves/Wavefield.hh"
#include "gz/waves/WavefieldSampler.hh"

#include <omp.h>
#include <random>

namespace gz
{
namespace waves
{

//////////////////////////////////////////////////
// Utilities
void DebugPrint(const geom::Triangle& triangle)
{
  gzmsg << "Vertex[0]:   " << triangle[0] << "\n";
  gzmsg << "Vertex[1]:   " << triangle[1] << "\n";
  gzmsg << "Vertex[2]:   " << triangle[2] << "\n";
  gzmsg << "Normal:      " << Geometry::Normal(triangle) << "\n";
}

//////////////////////////////////////////////////
// Physics

geom::Point3 Physics::CenterOfForce(
  double _fA, double _fB,
  const geom::Point3& _A,
  const geom::Point3& _B
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

//////////////////////////////////////////////////
double Physics::DeepWaterDispersionToOmega(double _wavenumber)
{
  const double g = std::fabs(PhysicalConstants::Gravity());
  return std::sqrt(g * _wavenumber);
}

//////////////////////////////////////////////////
double Physics::DeepWaterDispersionToWavenumber(double _omega)
{
  const double g = std::fabs(PhysicalConstants::Gravity());
  return _omega * _omega / g;
}

//////////////////////////////////////////////////
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

//////////////////////////////////////////////////
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

//////////////////////////////////////////////////
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

//////////////////////////////////////////////////
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

  // Initialise to the base midpoint (correct force calcuation for
  // triangles with a horizontal base)
  double fU = 0, fL = 0;
  geom::Point3 CpU = B;
  geom::Point3 CpL = B;

  // Upper triangle H > M
  if (_H.z() >= _M.z())
  {
    // Center of pressure
    double z0 = _depthC - (_H.z() - _C.z());
    CpU = CenterOfPressureApexUp(z0, _H, _M, B);

    // Force at centroid
    geom::Point3 CU = Geometry::TriangleCentroid(_H, _M, D);
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
    geom::Point3 CL = Geometry::TriangleCentroid(_L, _M, D);
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

//////////////////////////////////////////////////
void Physics::BuoyancyForceAtCentroid(
  const WavefieldSampler& _wavefieldSampler,
  const geom::Triangle& _triangle,
  geom::Point3& _center,
  geom::Vector3& _force
)
{
  // Physical constants
  double density = PhysicalConstants::WaterDensity();   // kg m^-3
  double gravity = PhysicalConstants::Gravity();        // m s^-1

  // Calculate the triangles centroid
  _center = Geometry::TriangleCentroid(_triangle);

  // Calculate the depth
  double h = _wavefieldSampler.ComputeDepth(_center);

  // Calculate the force
  geom::Vector3 normal = Geometry::Normal(_triangle);
  double area = Geometry::TriangleArea(_triangle);
  _force = normal * (density * gravity * area * h);
}

//////////////////////////////////////////////////
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
  geom::Point3 C = Geometry::TriangleCentroid(_triangle);

  // Calculate the depth
  double depthC = _wavefieldSampler.ComputeDepth(C);

  geom::Vector3 normal = Geometry::Normal(_triangle);

  // Calculate buoyancy
  BuoyancyForceAtCenterOfPressure(depthC, C, H, M, L, normal, _center, _force);
}

//////////////////////////////////////////////////
std::array<double, 3> Physics::ComputeHeightMap(
  const WavefieldSampler& _wavefieldSampler,
  const geom::Triangle& _triangle
)
{
  // Heightmap for the triangle vertices
  geom::Direction3 direction(0, 0, -1);
  std::array<double, 3> heightMap;

  // Calculate the height above the surface (-depth)
  for (Index i=0; i < 3; ++i)
  {
    geom::Point3 vertex = _triangle[i];
    heightMap[i] = -_wavefieldSampler.ComputeDepth(vertex);
  }
  return heightMap;
}

//////////////////////////////////////////////////
//////////////////////////////////////////////////
class HydrodynamicsParametersPrivate
{
 public:
  HydrodynamicsParametersPrivate() :
    dampingOn(true),
    cDampU1(1.0e-6), cDampU2(1.0e-6),
    cDampV1(1.0e-3), cDampV2(1.0e-3),
    cDampW1(1.0e-3), cDampW2(1.0e-3),
    cDampP1(5.0e-3), cDampP2(5.0e-3),
    cDampQ1(5.0e-3), cDampQ2(5.0e-3),
    cDampN1(5.0e-4), cDampN2(5.0e-4),
    viscousDragOn(true),
    pressureDragOn(true),
    cPDrag1(1.0E+2),
    cPDrag2(1.0E+2),
    fPDrag(0.4),
    cSDrag1(1.0E+2),
    cSDrag2(1.0E+2),
    fSDrag(0.4),
    vRDrag(1.0),
    foilLiftOn(true),
    cLift1(1.0),
    alphaStall(0.26), // ~15° in radians (tune per hull)
    cLMax(1.62193)
  {
  }

  // Linear and rotational damping
  bool dampingOn;

  // Per-DOF Fossen damping: surge(U) sway(V) heave(W) roll(P) pitch(Q) yaw(N)
  double cDampU1, cDampU2;  // surge
  double cDampV1, cDampV2;  // sway   (high: wide hull)
  double cDampW1, cDampW2;  // heave
  double cDampP1, cDampP2;  // roll   (Seakeeper approximation range)
  double cDampQ1, cDampQ2;  // pitch
  double cDampN1, cDampN2;  // yaw

  /// Viscous drag
  bool viscousDragOn;

  /// Pressure drag
  bool pressureDragOn;

  /// Positive pressure (lift)
  double cPDrag1;
  double cPDrag2;
  double fPDrag;
  /// Negative pressure (suction)
  double cSDrag1;
  double cSDrag2;
  double fSDrag;

  // Reference speed for the pressure drag calculation
  double vRDrag;

  static constexpr double defaultDampingDistMin = 1.0E-8;
  static constexpr double defaultDampingDistMax = 1.0E-3;
  static constexpr double defaultCPDragDistMin = 20.0;
  static constexpr double defaultCPDragDistMax = 400.0;
  static constexpr double defaultCSDragMin = 50.0;
  static constexpr double defaultCSDragMax = 800.0;
  static constexpr double defaultFPDragMin = 0.2;
  static constexpr double defaultFPDragMax = 0.7;
  static constexpr double defaultFSDragMin = 0.2;
  static constexpr double defaultFSDragMax = 0.7;
  static constexpr double defaultVRDragMin = 0.5;
  static constexpr double defaultVRDragMax = 2.0;

  bool   foilLiftOn;
  double cLift1;      // Cl scale factor (tune per hull)
  double alphaStall; // Stall angle of attack (rad)
  double cLMax;      // Maximum lift coefficient (tune per hull)

  WaterCurrentGrid water_current_grid_;
};

//////////////////////////////////////////////////
HydrodynamicsParameters::~HydrodynamicsParameters()
{
}

//////////////////////////////////////////////////
HydrodynamicsParameters::HydrodynamicsParameters() :
  data(new HydrodynamicsParametersPrivate())
{
}

//////////////////////////////////////////////////
bool HydrodynamicsParameters::DampingOn() const
{
  return this->data->dampingOn;
}

//////////////////////////////////////////////////
bool HydrodynamicsParameters::ViscousDragOn() const
{
  return this->data->viscousDragOn;
}

//////////////////////////////////////////////////
bool HydrodynamicsParameters::PressureDragOn() const
{
  return this->data->pressureDragOn;
}

//////////////////////////////////////////////////
double HydrodynamicsParameters::CDampU1() const { return this->data->cDampU1; }
double HydrodynamicsParameters::CDampU2() const { return this->data->cDampU2; }
double HydrodynamicsParameters::CDampV1() const { return this->data->cDampV1; }
double HydrodynamicsParameters::CDampV2() const { return this->data->cDampV2; }
double HydrodynamicsParameters::CDampW1() const { return this->data->cDampW1; }
double HydrodynamicsParameters::CDampW2() const { return this->data->cDampW2; }
double HydrodynamicsParameters::CDampP1() const { return this->data->cDampP1; }
double HydrodynamicsParameters::CDampP2() const { return this->data->cDampP2; }
double HydrodynamicsParameters::CDampQ1() const { return this->data->cDampQ1; }
double HydrodynamicsParameters::CDampQ2() const { return this->data->cDampQ2; }
double HydrodynamicsParameters::CDampN1() const { return this->data->cDampN1; }
double HydrodynamicsParameters::CDampN2() const { return this->data->cDampN2; }

//////////////////////////////////////////////////
double HydrodynamicsParameters::CPDrag1() const
{
  return this->data->cPDrag1;
}

//////////////////////////////////////////////////
double HydrodynamicsParameters::CPDrag2() const
{
  return this->data->cPDrag2;
}

//////////////////////////////////////////////////
double HydrodynamicsParameters::FPDrag() const
{
  return this->data->fPDrag;
}

//////////////////////////////////////////////////
double HydrodynamicsParameters::CSDrag1() const
{
  return this->data->cSDrag1;
}

//////////////////////////////////////////////////
double HydrodynamicsParameters::CSDrag2() const
{
  return this->data->cSDrag2;
}

//////////////////////////////////////////////////
double HydrodynamicsParameters::FSDrag() const
{
  return this->data->fSDrag;
}

//////////////////////////////////////////////////
double HydrodynamicsParameters::VRDrag() const
{
  return this->data->vRDrag;
}

//////////////////////////////////////////////////
bool HydrodynamicsParameters::FoilLiftOn() const
{
  return this->data->foilLiftOn;
}

//////////////////////////////////////////////////
double HydrodynamicsParameters::CLift1() const
{
  return this->data->cLift1;
}

//////////////////////////////////////////////////
double HydrodynamicsParameters::CLMax() const
{
  return this->data->cLMax;
}


//////////////////////////////////////////////////
double HydrodynamicsParameters::AlphaStall() const
{
  return this->data->alphaStall;
}

//////////////////////////////////////////////////
const WaterCurrentGrid& HydrodynamicsParameters::GetWaterCurrentGrid() const
{
  return this->data->water_current_grid_;
}

//////////////////////////////////////////////////
void HydrodynamicsParameters::SetFromMsg(const gz::msgs::Param_V& _msg)
{
  this->data->dampingOn      = Utilities::MsgParamBool(
      _msg,  "damping_on",       this->data->dampingOn);
  this->data->viscousDragOn  = Utilities::MsgParamBool(
      _msg,  "viscous_drag_on",  this->data->viscousDragOn);
  this->data->pressureDragOn = Utilities::MsgParamBool(
      _msg,  "pressure_drag_on", this->data->pressureDragOn);

  this->data->cDampU1 = Utilities::MsgParamDouble(_msg, "cDampU1", this->data->cDampU1);
  this->data->cDampU2 = Utilities::MsgParamDouble(_msg, "cDampU2", this->data->cDampU2);
  this->data->cDampV1 = Utilities::MsgParamDouble(_msg, "cDampV1", this->data->cDampV1);
  this->data->cDampV2 = Utilities::MsgParamDouble(_msg, "cDampV2", this->data->cDampV2);
  this->data->cDampW1 = Utilities::MsgParamDouble(_msg, "cDampW1", this->data->cDampW1);
  this->data->cDampW2 = Utilities::MsgParamDouble(_msg, "cDampW2", this->data->cDampW2);
  this->data->cDampP1 = Utilities::MsgParamDouble(_msg, "cDampP1", this->data->cDampP1);
  this->data->cDampP2 = Utilities::MsgParamDouble(_msg, "cDampP2", this->data->cDampP2);
  this->data->cDampQ1 = Utilities::MsgParamDouble(_msg, "cDampQ1", this->data->cDampQ1);
  this->data->cDampQ2 = Utilities::MsgParamDouble(_msg, "cDampQ2", this->data->cDampQ2);
  this->data->cDampN1 = Utilities::MsgParamDouble(_msg, "cDampN1", this->data->cDampN1);
  this->data->cDampN2 = Utilities::MsgParamDouble(_msg, "cDampN2", this->data->cDampN2);
  this->data->cPDrag1 = Utilities::MsgParamDouble(
      _msg, "cPDrag1",  this->data->cPDrag1);
  this->data->cPDrag2 = Utilities::MsgParamDouble(
      _msg, "cPDrag2",  this->data->cPDrag2);
  this->data->fPDrag  = Utilities::MsgParamDouble(
      _msg, "fPDrag",   this->data->fPDrag);
  this->data->cSDrag1 = Utilities::MsgParamDouble(
      _msg, "cSDrag1",  this->data->cSDrag1);
  this->data->cSDrag2 = Utilities::MsgParamDouble(
      _msg, "cSDrag2",  this->data->cSDrag2);
  this->data->fSDrag  = Utilities::MsgParamDouble(
      _msg, "fSDrag",   this->data->fSDrag);
  this->data->vRDrag  = Utilities::MsgParamDouble(
      _msg, "vRDrag",   this->data->vRDrag);

  this->data->foilLiftOn = Utilities::MsgParamBool(
      _msg, "foil_lift_on", this->data->foilLiftOn);
  this->data->cLift1 = Utilities::MsgParamDouble(
      _msg, "cLift1",     this->data->cLift1);
  this->data->alphaStall = Utilities::MsgParamDouble(
      _msg, "alphaStall", this->data->alphaStall);
  // cLMax default recomputed from (possibly updated) cLift1 and alphaStall
  this->data->cLMax = Utilities::MsgParamDouble(
      _msg, "cLMax",
      this->data->cLift1 * 2.0 * M_PI * std::sin(this->data->alphaStall));
}

//////////////////////////////////////////////////
void HydrodynamicsParameters::SetFromSDF(sdf::Element& _sdf)
{
  std::string bin_path;
  if (_sdf.HasElement("water_current_grid"))
    bin_path = _sdf.Get<std::string>("water_current_grid");

  if (!bin_path.empty())
    this->data->water_current_grid_.LoadFromFile(bin_path);

  this->data->foilLiftOn = Utilities::SdfParamBool(_sdf, "foil_lift_on", this->data->foilLiftOn);
  this->data->cLift1 = Utilities::SdfParamDouble(_sdf, "cLift1", this->data->cLift1);
  this->data->alphaStall = Utilities::SdfParamDouble(_sdf, "alphaStall", this->data->alphaStall);
  this->data->cLMax = Utilities::SdfParamDouble(_sdf, "cLMax",
    this->data->cLift1 * 2.0 * M_PI * std::sin(this->data->alphaStall));

  this->data->dampingOn = Utilities::SdfParamBool(
    _sdf, "damping_on", this->data->dampingOn);
  this->data->viscousDragOn = Utilities::SdfParamBool(
    _sdf, "viscous_drag_on", this->data->viscousDragOn);
  this->data->pressureDragOn = Utilities::SdfParamBool(
    _sdf, "pressure_drag_on", this->data->pressureDragOn);

  if (_sdf.HasElement("randomize"))
  {
    this->SetRandomFromSDF(_sdf);
    return;
  }

  this->data->cDampU1 = Utilities::SdfParamDouble(_sdf, "cDampU1", this->data->cDampU1);
  this->data->cDampU2 = Utilities::SdfParamDouble(_sdf, "cDampU2", this->data->cDampU2);
  this->data->cDampV1 = Utilities::SdfParamDouble(_sdf, "cDampV1", this->data->cDampV1);
  this->data->cDampV2 = Utilities::SdfParamDouble(_sdf, "cDampV2", this->data->cDampV2);
  this->data->cDampW1 = Utilities::SdfParamDouble(_sdf, "cDampW1", this->data->cDampW1);
  this->data->cDampW2 = Utilities::SdfParamDouble(_sdf, "cDampW2", this->data->cDampW2);
  this->data->cDampP1 = Utilities::SdfParamDouble(_sdf, "cDampP1", this->data->cDampP1);
  this->data->cDampP2 = Utilities::SdfParamDouble(_sdf, "cDampP2", this->data->cDampP2);
  this->data->cDampQ1 = Utilities::SdfParamDouble(_sdf, "cDampQ1", this->data->cDampQ1);
  this->data->cDampQ2 = Utilities::SdfParamDouble(_sdf, "cDampQ2", this->data->cDampQ2);
  this->data->cDampN1 = Utilities::SdfParamDouble(_sdf, "cDampN1", this->data->cDampN1);
  this->data->cDampN2 = Utilities::SdfParamDouble(_sdf, "cDampN2", this->data->cDampN2);
  this->data->cPDrag1 = Utilities::SdfParamDouble(
    _sdf, "cPDrag1", this->data->cPDrag1);
  this->data->cPDrag2 = Utilities::SdfParamDouble(
    _sdf, "cPDrag2", this->data->cPDrag2);
  this->data->fPDrag = Utilities::SdfParamDouble(
    _sdf, "fPDrag", this->data->fPDrag);
  this->data->cSDrag1 = Utilities::SdfParamDouble(
    _sdf, "cSDrag1", this->data->cSDrag1);
  this->data->cSDrag2 = Utilities::SdfParamDouble(
    _sdf, "cSDrag2", this->data->cSDrag2);
  this->data->fSDrag = Utilities::SdfParamDouble(
    _sdf, "fSDrag", this->data->fSDrag);
  this->data->vRDrag = Utilities::SdfParamDouble(
    _sdf, "vRDrag", this->data->vRDrag);
}

void HydrodynamicsParameters::SetRandomFromSDF(sdf::Element& _sdf) {
  // Seed the random number engine using the current time
  unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
  std::mt19937 engine(seed); // Using the Mersenne Twister 32-bit engine

  auto dampDist = [&](const char* minKey, const char* maxKey,
                      double defMin, double defMax) {
    return std::uniform_real_distribution<double>(
      Utilities::SdfParamDouble(_sdf, minKey, defMin),
      Utilities::SdfParamDouble(_sdf, maxKey, defMax));
  };
  auto dU = dampDist("dampUMin", "dampUMax", 1.0e-7, 1.0e-5);
  auto dV = dampDist("dampVMin", "dampVMax", 1.0e-4, 1.0e-2);
  auto dW = dampDist("dampWMin", "dampWMax", 1.0e-4, 1.0e-2);
  auto dP = dampDist("dampPMin", "dampPMax", 1.0e-3, 5.0e-2);
  auto dQ = dampDist("dampQMin", "dampQMax", 1.0e-3, 5.0e-2);
  auto dN = dampDist("dampNMin", "dampNMax", 1.0e-4, 1.0e-3);
  this->data->cDampU1 = dU(engine); this->data->cDampU2 = dU(engine);
  this->data->cDampV1 = dV(engine); this->data->cDampV2 = dV(engine);
  this->data->cDampW1 = dW(engine); this->data->cDampW2 = dW(engine);
  this->data->cDampP1 = dP(engine); this->data->cDampP2 = dP(engine);
  this->data->cDampQ1 = dQ(engine); this->data->cDampQ2 = dQ(engine);
  this->data->cDampN1 = dN(engine); this->data->cDampN2 = dN(engine);

  std::uniform_real_distribution<double> cPDragDist(
    Utilities::SdfParamDouble(_sdf, "cPDragDistMin", this->data->defaultCPDragDistMin),
    Utilities::SdfParamDouble(_sdf, "cPDragDistMax", this->data->defaultCPDragDistMax));
  this->data->cPDrag1 = cPDragDist(engine);
  this->data->cPDrag2 = cPDragDist(engine);

  std::uniform_real_distribution<double> fPDragDist(
    Utilities::SdfParamDouble(_sdf, "fPDragDistMin", this->data->defaultFPDragMin),
    Utilities::SdfParamDouble(_sdf, "fPDragDistMax", this->data->defaultFPDragMax));
  this->data->fPDrag = fPDragDist(engine);

  std::uniform_real_distribution<double> cSDragDist(
    Utilities::SdfParamDouble(_sdf, "cSDragDistMin", this->data->defaultCSDragMin),
    Utilities::SdfParamDouble(_sdf, "cSDragDistMax", this->data->defaultCSDragMax));
  this->data->cSDrag1 = cSDragDist(engine);
  this->data->cSDrag2 = cSDragDist(engine);

  std::uniform_real_distribution<double> fSDragDist(
    Utilities::SdfParamDouble(_sdf, "fSDragDistMin", this->data->defaultFSDragMin),
    Utilities::SdfParamDouble(_sdf, "fSDragDistMax", this->data->defaultFSDragMax));
  this->data->fSDrag = fSDragDist(engine);

  std::uniform_real_distribution<double> vRDragDist(
    Utilities::SdfParamDouble(_sdf, "vRDragDistMin", this->data->defaultVRDragMin),
    Utilities::SdfParamDouble(_sdf, "vRDragDistMax", this->data->defaultVRDragMax));
  this->data->vRDrag = vRDragDist(engine);
}

//////////////////////////////////////////////////
void HydrodynamicsParameters::DebugPrint() const
{
  gzmsg << "damping_on:       " << this->data->dampingOn << "\n";
  gzmsg << "viscous_drag_on:  " << this->data->viscousDragOn << "\n";
  gzmsg << "pressure_drag_on: " << this->data->pressureDragOn << "\n";
  gzmsg << "cDampU1/U2:       " << this->data->cDampU1 << " / " << this->data->cDampU2 << "\n";
  gzmsg << "cDampV1/V2:       " << this->data->cDampV1 << " / " << this->data->cDampV2 << "\n";
  gzmsg << "cDampW1/W2:       " << this->data->cDampW1 << " / " << this->data->cDampW2 << "\n";
  gzmsg << "cDampP1/P2:       " << this->data->cDampP1 << " / " << this->data->cDampP2 << "\n";
  gzmsg << "cDampQ1/Q2:       " << this->data->cDampQ1 << " / " << this->data->cDampQ2 << "\n";
  gzmsg << "cDampN1/N2:       " << this->data->cDampN1 << " / " << this->data->cDampN2 << "\n";
  gzmsg << "cPDrag1:          " << this->data->cPDrag1 << "\n";
  gzmsg << "cPDrag2:          " << this->data->cPDrag2 << "\n";
  gzmsg << "fPDrag:           " << this->data->fPDrag << "\n";
  gzmsg << "cSDrag1:          " << this->data->cSDrag1 << "\n";
  gzmsg << "cSDrag2:          " << this->data->cSDrag2 << "\n";
  gzmsg << "fSDrag:           " << this->data->fSDrag << "\n";
  gzmsg << "vRDrag:           " << this->data->vRDrag << "\n";
}

//////////////////////////////////////////////////
void DebugPrint(const TriangleProperties& props)
{
  gzmsg << "index:        " << props.index << "\n";
  gzmsg << "normal:       " << props.normal << "\n";
  gzmsg << "area:         " << props.area << "\n";
  gzmsg << "subArea:      " << props.subArea << "\n";
  gzmsg << "vh:           " << props.vh << "\n";
  gzmsg << "vm:           " << props.vm << "\n";
  gzmsg << "vl:           " << props.vl << "\n";
  gzmsg << "hh:           " << props.hh << "\n";
  gzmsg << "hm:           " << props.hm << "\n";
  gzmsg << "hl:           " << props.hl << "\n";
}

//////////////////////////////////////////////////
void DebugPrint(const SubmergedTriangleProperties& props)
{
  gzmsg << "index:        " << props.index << "\n";
  gzmsg << "normal:       " << props.normal << "\n";
  gzmsg << "centroid:     " << props.centroid << "\n";
  gzmsg << "xr:           " << props.xr << "\n";
  gzmsg << "area:         " << props.area << "\n";
  gzmsg << "vp:           " << props.vp << "\n";
  gzmsg << "up:           " << props.up << "\n";
  gzmsg << "cosTheta:     " << props.cosTheta << "\n";
  gzmsg << "vn:           " << props.vn << "\n";
  gzmsg << "vt:           " << props.vt << "\n";
  gzmsg << "ut:           " << props.ut << "\n";
  gzmsg << "uf:           " << props.uf << "\n";
  gzmsg << "vf:           " << props.vf << "\n";
}

//////////////////////////////////////////////////

class HydrodynamicsPrivate
{
 public:
  /// \brief The hydrodynamics parameters.
  std::shared_ptr<const HydrodynamicsParameters> params;

  /// \brief The mesh of the rigid body described by this model link.
  std::shared_ptr<const geom::Mesh> linkMesh;

  /// \brief The wavefield sampler for this rigid body (linkMesh).
  std::shared_ptr<const WavefieldSampler>  wavefieldSampler;

  /// \brief Pose of the centre of mass.
  gz::math::Pose3d pose;

  /// \brief Position of the centre of mass (CGAL types).
  geom::Point3 position;

  // \brief Linear velocity of the centre of mass.
  geom::Vector3 linVelocity;

  /// \brief Angular velocity of the centre of mass.
  geom::Vector3 angVelocity;

  /// \brief The calculated waterline length.
  double waterlineLength;

  /// \brief The calculated waterline beam (max width at waterline).
  double waterlineBeam;

  /// \brief The aspect ratio of the dynamic foil (for lift calculation).
  double dynamic_foil_ar;

  /// \brief The water current at the centre of mass.
  geom::Vector3 waterCurrentCoM;

  /// \brief The depth at each vertex point (indexed by vertex index).
  std::vector<double> depths;
  std::vector<geom::Triangle> submergedTriangles;
  std::vector<TriangleProperties> triangleProperties;
  std::vector<SubmergedTriangleProperties> submergedTriangleProperties;
  std::vector<geom::Line> waterline;

  /// \brief Per-thread scratch buffers for the parallel face loop in
  ///        UpdateSubmergedTriangles. Persisted across steps to avoid
  ///        malloc/free on every physics tick.
  std::vector<std::vector<geom::Triangle>>              tl_subTris;
  std::vector<std::vector<SubmergedTriangleProperties>> tl_subProps;
  std::vector<std::vector<geom::Line>>                  tl_waterlines;

  double area;

  double submergedArea;

  // Keep buoyance force and center of pressure for debugging...
  std::vector<geom::Vector3> fBuoyancy;
  std::vector<geom::Point3>  cBuoyancy;

  /// \brief The computed force
  geom::Vector3 force;

  /// \brief The computed torque
  geom::Vector3 torque;
};

//////////////////////////////////////////////////

Hydrodynamics::Hydrodynamics(
  std::shared_ptr<const HydrodynamicsParameters> _params,
  std::shared_ptr<const geom::Mesh> _linkMesh,
  std::shared_ptr<const WavefieldSampler> _wavefieldSampler
) : data(new HydrodynamicsPrivate())
{
  this->data->params = _params;
  this->data->linkMesh = _linkMesh;
  this->data->wavefieldSampler = _wavefieldSampler;
  this->data->position = geom::Origin();
  this->data->linVelocity = geom::NullVector();
  this->data->angVelocity = geom::NullVector();
  this->data->waterlineLength = 0.0;

  // Allocate the per-vertex depth array once; reused every physics step.
  this->data->depths.assign(geom::VertexCount(*_linkMesh), 0.0);
}

//////////////////////////////////////////////////
void Hydrodynamics::Update(
  std::shared_ptr<const WavefieldSampler> _wavefieldSampler,
  const gz::math::Pose3d& _pose,
  const geom::Vector3& _linVelocity,
  const geom::Vector3& _angVelocity,
  const std::chrono::_V2::steady_clock::duration& simTime
  )
{
  // Set rigid body props.
  this->data->wavefieldSampler = _wavefieldSampler;
  this->data->pose = _pose;
  this->data->position = ToPoint3(_pose.Pos());
  this->data->linVelocity = _linVelocity;
  this->data->angVelocity = _angVelocity;

  // Reset
  this->data->force = geom::NullVector();
  this->data->torque = geom::NullVector();

  // Update physics
  this->UpdateSubmergedTriangles();
  this->ComputeAreas();
  this->ComputeWaterlineLength();
  this->ComputeWaterlineBeam();
  this->ComputeDynamicFoilGeometry();
  this->SampleWaterCurrentCoM();

  this->ComputeAllSubmergedForces(simTime);

  if (this->data->params->DampingOn())
    this->ComputeDampingForce();
}

//////////////////////////////////////////////////
const geom::Vector3& Hydrodynamics::Force() const
{
  return this->data->force;
}

//////////////////////////////////////////////////
const geom::Vector3& Hydrodynamics::Torque() const
{
  return this->data->torque;
}

//////////////////////////////////////////////////
const std::vector<geom::Line>& Hydrodynamics::GetWaterline() const
{
  return this->data->waterline;
}

//////////////////////////////////////////////////
const std::vector<geom::Triangle>& Hydrodynamics::GetSubmergedTriangles() const
{
  return this->data->submergedTriangles;
}

const std::vector<TriangleProperties>& Hydrodynamics::GetTriangleProperties() const
{
  return this->data->triangleProperties;
}

const std::vector<SubmergedTriangleProperties>& Hydrodynamics::GetSubmergedTriangleProperties() const
{
  return this->data->submergedTriangleProperties;
}

const gz::waves::geom::Vector3 Hydrodynamics::GetWaterCurrentCoM() const
{
  return this->data->waterCurrentCoM;
}

//////////////////////////////////////////////////
void Hydrodynamics::UpdateSubmergedTriangles()
{
  // submergedTriangles/submergedTriangleProperties/waterline are rebuilt every
  // step; triangleProperties is resize()-only so it skips re-init when nFaces
  // is unchanged (hull mesh is fixed at runtime).
  this->data->submergedTriangles.clear();
  this->data->submergedTriangleProperties.clear();
  this->data->waterline.clear();

  auto& linkMesh = *this->data->linkMesh;
  auto& wavefieldSampler = *this->data->wavefieldSampler;

  // Compute depths — depth array allocated once in constructor.
  // Vertex/face indices are 0..N-1 (mesh topology never changes, no deletions).
  const int nVerts = static_cast<int>(geom::VertexCount(linkMesh));
  const int nTV = std::max(1, std::min(omp_get_max_threads(), nVerts / 32));
  #pragma omp parallel for schedule(static) num_threads(nTV)
  for (int i = 0; i < nVerts; ++i)
  {
    this->data->depths[i] =
        wavefieldSampler.ComputeDepth(geom::VertexPoint(linkMesh, i));
  }

  const int nFaces = static_cast<int>(geom::FaceCount(linkMesh));
  this->data->triangleProperties.resize(nFaces);

  // Thread-local output buffers — threads push_back independently, merged below.
  // Cap threads so that each thread gets at least 32 faces; with fewer faces
  // per thread the OpenMP barrier overhead exceeds the per-iteration work.
  const int nThreads = std::max(1, std::min(omp_get_max_threads(), nFaces / 32));

  // Resize persistent per-thread buffers only when thread count changes, then
  // clear each step. This avoids malloc/free on every physics tick.
  auto& tl_subTris    = this->data->tl_subTris;
  auto& tl_subProps   = this->data->tl_subProps;
  auto& tl_waterlines = this->data->tl_waterlines;
  if (static_cast<int>(tl_subTris.size()) != nThreads)
  {
    const int reservePerThread = (nFaces / nThreads) * 2 + 4;
    tl_subTris.resize(nThreads);
    tl_subProps.resize(nThreads);
    tl_waterlines.resize(nThreads);
    for (int t = 0; t < nThreads; ++t)
    {
      tl_subTris[t].reserve(reservePerThread);
      tl_subProps[t].reserve(reservePerThread);
      tl_waterlines[t].reserve(reservePerThread / 2);
    }
  }
  for (int t = 0; t < nThreads; ++t)
  {
    tl_subTris[t].clear();
    tl_subProps[t].clear();
    tl_waterlines[t].clear();
  }

  #pragma omp parallel for schedule(static) num_threads(nThreads)
  for (int i = 0; i < nFaces; ++i)
  {
    const int tid = omp_get_thread_num();
    geom::Triangle triangle = geom::FaceTriangle(linkMesh, i);

    TriangleProperties& triProps = this->data->triangleProperties[i];
    triProps.normal = Geometry::Normal(triangle);
    triProps.area   = Geometry::TriangleArea(triangle);

    // Note sign change for height.
    const auto fv = geom::FaceVertices(linkMesh, i);
    for (int j = 0; j < 3; ++j)
      triProps.heightMap[j] = -this->data->depths[fv[j]];

    this->PopulateSubmergedTriangle(
        triangle, triProps,
        tl_subTris[tid], tl_subProps[tid], tl_waterlines[tid]);
  }

  // Serial merge of per-thread results into shared data.
  for (int t = 0; t < nThreads; ++t)
  {
    this->data->submergedTriangles.insert(
        this->data->submergedTriangles.end(),
        std::make_move_iterator(tl_subTris[t].begin()),
        std::make_move_iterator(tl_subTris[t].end()));
    this->data->submergedTriangleProperties.insert(
        this->data->submergedTriangleProperties.end(),
        std::make_move_iterator(tl_subProps[t].begin()),
        std::make_move_iterator(tl_subProps[t].end()));
    this->data->waterline.insert(
        this->data->waterline.end(),
        std::make_move_iterator(tl_waterlines[t].begin()),
        std::make_move_iterator(tl_waterlines[t].end()));
  }
}

//////////////////////////////////////////////////
void Hydrodynamics::PopulateSubmergedTriangle(
  const geom::Triangle& _triangle,
  TriangleProperties& _triProps,
  std::vector<geom::Triangle>& _subTris,
  std::vector<SubmergedTriangleProperties>& _subProps,
  std::vector<geom::Line>& _waterlines)
{
  // Calculations
  const Index H = 0, M = 1, L = 2;
  std::array<Index, 3> idx = algorithm::sort_indexes(_triProps.heightMap);

  _triProps.hh = _triProps.heightMap[idx[H]];
  _triProps.hm = _triProps.heightMap[idx[M]];
  _triProps.hl = _triProps.heightMap[idx[L]];

  _triProps.vh = _triangle[idx[H]];
  _triProps.vm = _triangle[idx[M]];
  _triProps.vl = _triangle[idx[L]];

  if (_triProps.hh > 0)
  {
    if (_triProps.hm > 0)
    {
      if (_triProps.hl > 0)
      {
        // no-op
      } else {
        this->SplitPartiallySubmergedTriangle1(_triProps, _subTris, _subProps, _waterlines);
      }
    } else {
      this->SplitPartiallySubmergedTriangle2(_triProps, _subTris, _subProps, _waterlines);
    }
  } else {
    this->AddFullySubmergedTriangle(_triProps, _subTris, _subProps);
  }
}

//////////////////////////////////////////////////
void Hydrodynamics::SplitPartiallySubmergedTriangle1(
    TriangleProperties& _triProps,
    std::vector<geom::Triangle>& _subTris,
    std::vector<SubmergedTriangleProperties>& _subProps,
    std::vector<geom::Line>& _waterlines)
{
  geom::Vector3& n = _triProps.normal;
  geom::Point3& vh = _triProps.vh;
  geom::Point3& vm = _triProps.vm;
  geom::Point3& vl = _triProps.vl;
  double hh = _triProps.hh;
  double hm = _triProps.hm;
  double hl = _triProps.hl;

  double tm = -hl/(hm - hl);
  double th = -hl/(hh - hl);

  geom::Point3 vmi = vl + (vm - vl) * tm;
  geom::Point3 vhi = vl + (vh - vl) * th;

  // Create the new submerged triangle
  geom::Triangle tri0(vl, vmi, vhi);
  if (geom::Dot(n, Geometry::Normal(tri0)) < 0.0)
  {
    // Change orientation
    tri0 = geom::Triangle(vl, vhi, vmi);
  }
  _subTris.push_back(tri0);

  // Properties of the submerged tri0
  SubmergedTriangleProperties subTriProps0;
  subTriProps0.index = _triProps.index;
  subTriProps0.normal = Geometry::Normal(tri0);
  subTriProps0.centroid = Geometry::TriangleCentroid(tri0);
  subTriProps0.area = Geometry::TriangleArea(tri0);
  _subProps.push_back(subTriProps0);

  // Fraction of original triangle submerged.
  _triProps.subArea = subTriProps0.area;

  // Create a new line (for the water line)
  _waterlines.emplace_back(vmi, vhi);
}

//////////////////////////////////////////////////
void Hydrodynamics::SplitPartiallySubmergedTriangle2(
    TriangleProperties& _triProps,
    std::vector<geom::Triangle>& _subTris,
    std::vector<SubmergedTriangleProperties>& _subProps,
    std::vector<geom::Line>& _waterlines)
{
  geom::Vector3& n = _triProps.normal;
  geom::Point3& vh = _triProps.vh;
  geom::Point3& vm = _triProps.vm;
  geom::Point3& vl = _triProps.vl;
  double hh = _triProps.hh;
  double hm = _triProps.hm;
  double hl = _triProps.hl;

  double tm = -hm/(hh - hm);
  double tl = -hl/(hh - hl);

  geom::Point3 vmi =  vm + (vh - vm) * tm;
  geom::Point3 vli =  vl + (vh - vl) * tl;

  // Create the new submerged triangles
  geom::Triangle tri0(vm, vmi, vl);
  geom::Triangle tri1(vmi, vli, vl);

  if (geom::Dot(n, Geometry::Normal(tri0)) < 0.0)
  {
    tri0 = geom::Triangle(vmi, vm, vl);
  }
  if (geom::Dot(n, Geometry::Normal(tri1)) < 0.0)
  {
    tri1 = geom::Triangle(vli, vmi, vl);
  }

  _subTris.push_back(tri0);
  _subTris.push_back(tri1);

  // Properties of the submerged tri0
  SubmergedTriangleProperties subTriProps0;
  subTriProps0.index = _triProps.index;
  subTriProps0.normal = Geometry::Normal(tri0);
  subTriProps0.centroid = Geometry::TriangleCentroid(tri0);
  subTriProps0.area = Geometry::TriangleArea(tri0);
  _subProps.push_back(subTriProps0);

  // Properties of the submerged tri1
  SubmergedTriangleProperties subTriProps1;
  subTriProps1.index = _triProps.index;
  subTriProps1.normal = Geometry::Normal(tri1);
  subTriProps1.centroid = Geometry::TriangleCentroid(tri1);
  subTriProps1.area = Geometry::TriangleArea(tri1);
  _subProps.push_back(subTriProps1);

  // Fraction of original triangle submerged.
  _triProps.subArea = subTriProps0.area + subTriProps1.area;

  // Create a new line (for the water line)
  _waterlines.emplace_back(vmi, vli);
}

//////////////////////////////////////////////////
void Hydrodynamics::AddFullySubmergedTriangle(
    TriangleProperties& _triProps,
    std::vector<geom::Triangle>& _subTris,
    std::vector<SubmergedTriangleProperties>& _subProps)
{
  // Add the full triangle
  geom::Vector3& n = _triProps.normal;
  geom::Point3& vh = _triProps.vh;
  geom::Point3& vm = _triProps.vm;
  geom::Point3& vl = _triProps.vl;

  // Create the new submerged triangle
  geom::Triangle tri(vh, vm, vl);
  if (geom::Dot(n, Geometry::Normal(tri)) < 0.0)
  {
    tri = geom::Triangle(vm, vh, vl);
  }
  _subTris.push_back(tri);

  // Properties of the submerged triangle
  SubmergedTriangleProperties subTriProps;
  subTriProps.index = _triProps.index;
  subTriProps.normal = _triProps.normal;
  subTriProps.centroid = Geometry::TriangleCentroid(tri);
  subTriProps.area = _triProps.area;
  _subProps.push_back(subTriProps);

  // Fraction of original triangle submerged.
  _triProps.subArea = subTriProps.area;
}

//////////////////////////////////////////////////
void Hydrodynamics::ComputeAreas()
{
  double area = 0.0;
  for (auto&& props : this->data->triangleProperties)
  {
    area += props.area;
  }
  this->data->area = area;

  double subArea = 0.0;
  for (auto&& props : this->data->submergedTriangleProperties)
  {
    subArea += props.area;
  }
  this->data->submergedArea = subArea;
}

//////////////////////////////////////////////////
void Hydrodynamics::ComputeWaterlineLength()
{
  // Calculate the direction of the x-axis
  geom::Vector3 xaxis = ToVector3(this->data->pose.Rot().RotateVector(
    gz::math::Vector3d(1, 0, 0)));

  if (this->data->waterline.empty())
  {
    this->data->waterlineLength = 0.0;
    return;
  }

  // Exact LWL: longitudinal extent of all waterline endpoints.
  // Handles non-convex hulls correctly; no 0.5 approximation needed.
  double minProj = std::numeric_limits<double>::max();
  double maxProj = std::numeric_limits<double>::lowest();
  for (auto&& line : this->data->waterline)
  {
    double p0 = geom::ToDouble(
        geom::Dot(line.point() - geom::Origin(), xaxis));
    double p1 = p0 + geom::ToDouble(
        geom::Dot(line.to_vector(), xaxis));
    minProj = std::min(minProj, std::min(p0, p1));
    maxProj = std::max(maxProj, std::max(p0, p1));
  }
  this->data->waterlineLength = maxProj - minProj;
  // @DEBUG_INFO
  // gzmsg << "waterline length: " << this->data->waterlineLength << "\n";
}

void Hydrodynamics::ComputeWaterlineBeam()
{
    // Reuse the same waterline loop, project onto y-axis instead
    geom::Vector3 yaxis = ToVector3(this->data->pose.Rot().RotateVector(
        gz::math::Vector3d(0, 1, 0)));

    if (this->data->waterline.empty())
    {
        this->data->waterlineBeam = 0.0;
        return;
    }

    // Exact beam: lateral extent of all waterline endpoints.
    double minProj = std::numeric_limits<double>::max();
    double maxProj = std::numeric_limits<double>::lowest();
    for (auto&& line : this->data->waterline)
    {
        double p0 = geom::ToDouble(
            geom::Dot(line.point() - geom::Origin(), yaxis));
        double p1 = p0 + geom::ToDouble(
            geom::Dot(line.to_vector(), yaxis));
        minProj = std::min(minProj, std::min(p0, p1));
        maxProj = std::max(maxProj, std::max(p0, p1));
    }
    this->data->waterlineBeam = maxProj - minProj;
}

void Hydrodynamics::ComputeDynamicFoilGeometry()
{
  // ── identify bottom triangles ─────────────────────────────────────────
  // "Bottom" = outward normal has significant upward Z component
  // (hull bottom normals point downward into water, so outward = upward in world)
  double wetted_bottom_area = 0.0;
  for (auto& props : this->data->submergedTriangleProperties)
  {
    double nz = geom::ToDouble(props.normal.z());
    if (nz > this->data->params->BOTTOM_THRESHOLD)  // upward-facing = bottom surface
      wetted_bottom_area += props.area;
  }

  double span = this->data->waterlineBeam;    // from ComputeWaterlineBeam()

  this->data->dynamic_foil_ar =
    (span * span) / (wetted_bottom_area + 1e-9);

  // @DEBUG_INFO
  // gzmsg << "Submerged area: " << this->data->submergedArea << "\n";
  // gzmsg << "wetted_bottom_area: " << wetted_bottom_area << "\n";
  // gzmsg << "span: " << span << "\n";
  // gzmsg << "dynamic_foil_ar: " << this->data->dynamic_foil_ar << "\n";
}


//////////////////////////////////////////////////
// Compute the point velocity at a triangles centroid
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
  props.v_rel_mag = std::sqrt(geom::ToDouble(geom::SquaredLength(props.v_rel)));

  const double v_rel_dot_n = geom::ToDouble(
      geom::Dot(props.v_rel, props.normal));
  props.v_rel_n = props.normal * v_rel_dot_n;
  props.v_rel_t = props.v_rel - props.v_rel_n;

  const double v_rel_t_mag = std::sqrt(
      geom::ToDouble(geom::SquaredLength(props.v_rel_t)));
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

///////////////////////////////////////////////////
void Hydrodynamics::SampleWaterCurrentCoM() {
  this->data->waterCurrentCoM = this->data->params->GetWaterCurrentGrid().
    SampleAt(this->data->position.x(), this->data->position.y());
}


//////////////////////////////////////////////////
// Compute the Reynolds number
double Hydrodynamics::ComputeReynoldsNumber() const
{
  // fluid speed
  geom::Vector3 v_rel = this->data->linVelocity - this->data->waterCurrentCoM;
  double u = std::sqrt(geom::SquaredLength(v_rel));

  // characteristic length
  double L = this->data->waterlineLength;

  // Reynolds number
  double kv = PhysicalConstants::WaterKinematicViscosity();
  double Rn = u * L / kv;

  return Rn;
}

//////////////////////////////////////////////////
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

//////////////////////////////////////////////////
void Hydrodynamics::ComputeDampingForce()
{
  auto& params = *this->data->params;

  // Transform velocities to body frame for per-DOF Fossen damping matrix.
  // The gz::math quaternion is the world-to-body rotation stored in pose.Rot().
  gz::math::Quaterniond R     = this->data->pose.Rot();
  gz::math::Quaterniond R_inv = R.Inverse();

  gz::math::Vector3d linW = ToGz(geom::Vector3(this->data->linVelocity
                                  - this->data->waterCurrentCoM));
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

//////////////////////////////////////////////////
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

//////////////////////////////////////////////////
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
      geom::ToDouble(geom::SquaredLength(lift_dir)));
  if (ld_mag < 1e-9)
    return {geom::NullVector(), geom::NullVector()};

  lift_dir = lift_dir / ld_mag;
  const geom::Vector3 F_foil =
      lift_dir * (Cl * q_A) + (-props.up) * (Cdi * q_A);
  return {F_foil, geom::Cross(props.xr, F_foil)};
}

void Hydrodynamics::ComputeAllSubmergedForces(
    const std::chrono::_V2::steady_clock::duration& simTime)
{
  const int n = static_cast<int>(this->data->submergedTriangleProperties.size());
  if (n == 0) return;

  // Resize only — every element is overwritten in the loop below, so no init needed.
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

}  // namespace waves
}  // namespace gz
