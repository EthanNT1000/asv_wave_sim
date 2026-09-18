// Copyright (C) 2026  Rhys Mainwaring and contributors
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

/// \file geom_baseline.cc
/// \brief Prints golden values for tests/geom/test_geom_baseline.py.
///
/// Every quantity here goes through the geometry facade (geom::) and the
/// force pipeline of Physics.cc, so it pins the behaviour that the CGAL
/// removal phases must preserve:
///   * wave-surface height at sample points (regular, trochoid and FFT
///     wavefields) via Wavefield::Height,
///   * depth below the sampled water patch via WavefieldSampler::ComputeDepth,
///   * hydrodynamic force / torque, submerged area, displaced volume and
///     submerged-triangle count for a box hull at several poses and
///     velocities via Hydrodynamics::Update,
///   * first ray/mesh intersection via Geometry::MakeAABBTree/SearchMesh.
///
/// Usage: geom_baseline [--bench N]   (JSON on stdout)

#include <chrono>
#include <cmath>
#include <cstdio>
#include <cstring>
#include <iomanip>
#include <iostream>
#include <memory>
#include <random>
#include <sstream>
#include <string>
#include <vector>

#include <gz/common/MeshManager.hh>
#include <gz/math/Pose3.hh>
#include <gz/math/Vector3.hh>

#include "gz/waves/Convert.hh"
#include "gz/waves/Geometry.hh"
#include "gz/waves/Grid.hh"
#include "gz/waves/MeshTools.hh"
#include "gz/waves/PhysicalConstants.hh"
#include "gz/waves/Physics.hh"
#include "gz/waves/Wavefield.hh"
#include "gz/waves/WavefieldSampler.hh"
#include "gz/waves/WaveParameters.hh"
#include "gz/waves/geom/Geom.hh"

using namespace gz;
using namespace gz::waves;

namespace
{
std::string Num(double v)
{
  std::ostringstream os;
  os << std::setprecision(17) << v;
  std::string s = os.str();
  if (s == "-0") s = "0";
  return s;
}

std::string Vec(const geom::Vector3& v)
{
  return "[" + Num(v.x()) + ", " + Num(v.y()) + ", " + Num(v.z()) + "]";
}

std::string Pt(const geom::Point3& p)
{
  return "[" + Num(p.x()) + ", " + Num(p.y()) + ", " + Num(p.z()) + "]";
}

std::shared_ptr<geom::Mesh> BoxMesh(const std::string& name,
    const math::Vector3d& size)
{
  auto* mgr = common::MeshManager::Instance();
  if (!mgr->HasMesh(name))
    mgr->CreateBox(name, size, math::Vector2d(1, 1));
  auto mesh = std::make_shared<geom::Mesh>();
  MeshTools::MakeSurfaceMesh(*mgr->MeshByName(name), *mesh);
  return mesh;
}

std::shared_ptr<WaveParameters> MakeParams(const std::string& algorithm)
{
  auto params = std::make_shared<WaveParameters>();
  params->SetAlgorithm(algorithm);
  params->SetTileSize(64.0);
  params->SetCellCount(64);
  params->SetNumber(3);
  params->SetAngle(0.6);
  params->SetScale(1.2);
  params->SetSteepness(1.0);
  params->SetAmplitude(1.5);
  params->SetPeriod(6.0);
  params->SetDirection(math::Vector2d(1.0, 0.3));
  params->SetWindSpeedAndAngle(8.0, 0.4);
  return params;
}

// Sample points spread over the tile, some near cell edges/diagonals.
std::vector<geom::Point3> SamplePoints()
{
  std::vector<geom::Point3> pts;
  const double xs[] = {-30.0, -12.5, -3.0, 0.0, 0.5, 1.0, 7.25, 19.9, 31.9};
  const double ys[] = {-25.0, -1.0, 0.0, 2.5, 4.0, 13.7, 30.1};
  for (double x : xs)
    for (double y : ys)
      pts.emplace_back(x, y, -0.7);
  return pts;
}
}  // namespace

int main(int argc, char** argv)
{
  int bench = 0;
  for (int i = 1; i < argc; ++i)
  {
    if (std::strcmp(argv[i], "--bench") == 0 && i + 1 < argc)
      bench = std::atoi(argv[++i]);
  }

  std::ostringstream out;
  out << "{\n";

  // ---------------------------------------------------------------------
  // Wave heights
  // ---------------------------------------------------------------------
  out << "  \"wave_height\": {\n";
  // "trochoid" is excluded: TrochoidIrregularWaveSimulation::ElevationAt
  // assigns a zero-column array into an Eigen::Ref and trips Eigen's resize
  // assertion in assert-enabled builds (pre-existing, unrelated to geom::).
  const char* algorithms[] = {"sinusoid", "fft"};
  const double times[] = {0.0, 1.7, 5.3};
  const auto pts = SamplePoints();
  for (int a = 0; a < 2; ++a)
  {
    auto wavefield = std::make_shared<Wavefield>("geom_baseline");
    wavefield->SetParameters(MakeParams(algorithms[a]));
    out << "    \"" << algorithms[a] << "\": {\n";
    for (int t = 0; t < 3; ++t)
    {
      wavefield->Update(times[t]);
      out << "      \"t=" << times[t] << "\": [";
      for (size_t i = 0; i < pts.size(); ++i)
      {
        double h = 0.0;
        bool ok = wavefield->Height(
            Eigen::Vector3d(pts[i].x(), pts[i].y(), pts[i].z()), h);
        out << (i ? ", " : "") << (ok ? Num(h) : std::string("null"));
      }
      out << "]" << (t < 2 ? "," : "") << "\n";
    }
    out << "    }" << (a < 1 ? "," : "") << "\n";
  }
  out << "  },\n";

  // ---------------------------------------------------------------------
  // Hydrodynamics on a box hull in an FFT sea and in still water
  // ---------------------------------------------------------------------
  auto wavefield = std::make_shared<Wavefield>("geom_baseline");
  wavefield->SetParameters(MakeParams("fft"));
  wavefield->Update(1.7);

  auto still = std::make_shared<Wavefield>("geom_baseline");
  auto stillParams = MakeParams("sinusoid");
  stillParams->SetAmplitude(0.0);
  still->SetParameters(stillParams);
  still->Update(0.0);

  struct Case
  {
    const char* name;
    math::Vector3d size;
    math::Pose3d pose;
    math::Vector3d linVel;
    math::Vector3d angVel;
    bool waves;
  };
  const Case cases[] = {
    {"still_box_1x1x1_rest",  {1, 1, 1},  {0, 0, 0, 0, 0, 0}, {0, 0, 0}, {0, 0, 0}, false},
    {"still_box_10x4x2_rest", {10, 4, 2}, {0, 0, 0, 0, 0, 0}, {0, 0, 0}, {0, 0, 0}, false},
    {"still_box_10x4x2_heel", {10, 4, 2}, {0.3, -0.2, 0.15, 0.12, -0.05, 0.4}, {0, 0, 0}, {0, 0, 0}, false},
    {"still_box_10x4x2_moving", {10, 4, 2}, {0.3, -0.2, 0.15, 0.12, -0.05, 0.4}, {2.0, 0.5, -0.3}, {0.1, -0.2, 0.05}, false},
    {"waves_box_10x4x2_rest", {10, 4, 2}, {1.0, 2.0, -0.2, 0.05, 0.1, 0.7}, {0, 0, 0}, {0, 0, 0}, true},
    {"waves_box_10x4x2_moving", {10, 4, 2}, {1.0, 2.0, -0.2, 0.05, 0.1, 0.7}, {3.0, -1.0, 0.2}, {-0.15, 0.1, 0.3}, true},
  };

  const double rho = PhysicalConstants::WaterDensity();
  const double g = std::fabs(PhysicalConstants::Gravity());

  out << "  \"hydrodynamics\": {\n";
  int c = 0;
  const int nCases = sizeof(cases) / sizeof(cases[0]);
  for (const auto& cs : cases)
  {
    std::string meshName = "geom_baseline_box_" + std::to_string(c);
    auto initMesh = BoxMesh(meshName, cs.size);
    auto mesh = std::make_shared<geom::Mesh>(*initMesh);

    // Apply the pose to the mesh (as Hydrodynamics.cc::ApplyPose does).
    const Index nV = geom::VertexCount(*initMesh);
    for (Index i = 0; i < nV; ++i)
    {
      const geom::Point3& p0 = geom::VertexPoint(*initMesh, i);
      math::Vector3d p1 = cs.pose.Rot().RotateVector(
          {p0.x(), p0.y(), p0.z()}) + cs.pose.Pos();
      geom::SetVertexPoint(*mesh, i, geom::Point3(p1.X(), p1.Y(), p1.Z()));
    }

    const double patchSize = 2.0 * std::max({cs.size.X(), cs.size.Y(),
        cs.size.Z(), 1.0});
    auto patch = std::make_shared<Grid>(
        std::array<double, 2>{patchSize, patchSize},
        std::array<Index, 2>{4, 4});
    auto sampler = std::make_shared<WavefieldSampler>(
        cs.waves ? wavefield : still, patch);
    sampler->ApplyPose(cs.pose);
    sampler->UpdatePatch();

    auto params = std::make_shared<HydrodynamicsParameters>();
    Hydrodynamics hydro(params, mesh, sampler);
    hydro.Update(sampler, cs.pose, ToVector3(cs.linVel), ToVector3(cs.angVel),
        std::chrono::duration_cast<std::chrono::steady_clock::duration>(
            std::chrono::duration<double>(1.7)));

    double subArea = 0.0;
    for (const auto& props : hydro.GetSubmergedTriangleProperties())
      subArea += props.area;

    // Depth below the patch at the hull vertices (exercises ComputeDepth).
    double depthSum = 0.0;
    for (Index i = 0; i < nV; ++i)
      depthSum += sampler->ComputeDepth(geom::VertexPoint(*mesh, i));

    out << "    \"" << cs.name << "\": {\n";
    out << "      \"force\": " << Vec(hydro.Force()) << ",\n";
    out << "      \"torque\": " << Vec(hydro.Torque()) << ",\n";
    out << "      \"submerged_area\": " << Num(subArea) << ",\n";
    out << "      \"submerged_triangles\": "
        << hydro.GetSubmergedTriangles().size() << ",\n";
    out << "      \"waterline_segments\": " << hydro.GetWaterline().size()
        << ",\n";
    out << "      \"vertex_depth_sum\": " << Num(depthSum) << ",\n";
    // Displaced volume implied by the vertical hydrostatic force at rest.
    out << "      \"displaced_volume_from_fz\": "
        << Num(hydro.Force().z() / (rho * g)) << "\n";
    out << "    }" << (++c < nCases ? "," : "") << "\n";
  }
  out << "  },\n";

  // ---------------------------------------------------------------------
  // Ray / mesh first intersection (AABB path)
  // ---------------------------------------------------------------------
  {
    auto mesh = BoxMesh("geom_baseline_box_ray", {10, 4, 2});
    auto tree = Geometry::MakeAABBTree(*mesh);
    struct RayCase { geom::Point3 o; geom::Direction3 d; };
    const RayCase rays[] = {
      {{0, 0, 0}, {0, 0, 1}},
      {{0, 0, 0}, {0, 0, -1}},
      {{1.5, -0.5, 0.3}, {1, 0, 0}},
      {{1.5, -0.5, 0.3}, {0.2, 0.7, -0.1}},
      {{20, 0, 0}, {-1, 0, 0}},
      {{0, 0, 5}, {0, 0, 1}},
    };
    out << "  \"ray_mesh\": [\n";
    int r = 0;
    for (const auto& rc : rays)
    {
      geom::Point3 hit = geom::Origin();
      bool ok = Geometry::SearchMesh(*tree, rc.o, rc.d, hit);
      out << "    " << (ok ? Pt(hit) : std::string("null"))
          << (++r < 6 ? "," : "") << "\n";
    }
    out << "  ]";
  }

  // ---------------------------------------------------------------------
  // Optional benchmark of the ray/mesh path and the depth (hot) path.
  // ---------------------------------------------------------------------
  if (bench > 0)
  {
    // A closed 4 x 2 x 1 m box tessellated 8 x 8 per face (768 triangles,
    // comparable to the 406/812-face WAM-V collision meshes), built directly
    // so it contains no degenerate triangles.
    auto sphere = std::make_shared<geom::Mesh>();
    {
      const int n = 8;
      const double hx = 2.0, hy = 1.0, hz = 0.5;
      auto addFace = [&](geom::Vector3 u, geom::Vector3 v, geom::Vector3 o)
      {
        std::vector<Index> idx;
        for (int i = 0; i <= n; ++i)
          for (int j = 0; j <= n; ++j)
          {
            double a = -1.0 + 2.0 * i / n, b = -1.0 + 2.0 * j / n;
            geom::Vector3 pv = o + u * a + v * b;
            idx.push_back(geom::AddVertex(*sphere,
                geom::Point3(pv.x(), pv.y(), pv.z())));
          }
        for (int i = 0; i < n; ++i)
          for (int j = 0; j < n; ++j)
          {
            Index i00 = idx[i * (n + 1) + j], i10 = idx[(i + 1) * (n + 1) + j];
            Index i01 = idx[i * (n + 1) + j + 1];
            Index i11 = idx[(i + 1) * (n + 1) + j + 1];
            geom::AddFace(*sphere, i00, i10, i11);
            geom::AddFace(*sphere, i00, i11, i01);
          }
      };
      addFace({hx, 0, 0}, {0, hy, 0}, {0, 0, hz});    // +z
      addFace({0, hy, 0}, {hx, 0, 0}, {0, 0, -hz});   // -z
      addFace({0, hy, 0}, {0, 0, hz}, {hx, 0, 0});    // +x
      addFace({0, 0, hz}, {0, hy, 0}, {-hx, 0, 0});   // -x
      addFace({0, 0, hz}, {hx, 0, 0}, {0, hy, 0});    // +y
      addFace({hx, 0, 0}, {0, 0, hz}, {0, -hy, 0});   // -y
    }

    std::mt19937 rng(12345);
    std::uniform_real_distribution<double> u(-1.0, 1.0);
    std::vector<geom::Point3> origins;
    std::vector<geom::Direction3> dirs;
    for (int i = 0; i < bench; ++i)
    {
      origins.emplace_back(0.5 * u(rng), 0.5 * u(rng), 0.5 * u(rng));
      dirs.emplace_back(u(rng), u(rng), u(rng));
    }

    auto t0 = std::chrono::steady_clock::now();
    auto tree = Geometry::MakeAABBTree(*sphere);
    auto t1 = std::chrono::steady_clock::now();
    int hits = 0;
    geom::Point3 hit = geom::Origin();
    for (int i = 0; i < bench; ++i)
      hits += Geometry::SearchMesh(*tree, origins[i], dirs[i], hit) ? 1 : 0;
    auto t2 = std::chrono::steady_clock::now();

    // Hot path: depth of random points below a wave patch.
    auto patch = std::make_shared<Grid>(
        std::array<double, 2>{20.0, 20.0}, std::array<Index, 2>{4, 4});
    auto sampler = std::make_shared<WavefieldSampler>(wavefield, patch);
    sampler->ApplyPose(math::Pose3d());
    sampler->UpdatePatch();
    std::vector<geom::Point3> queries;
    for (int i = 0; i < bench; ++i)
      queries.emplace_back(9.0 * u(rng), 9.0 * u(rng), u(rng));
    auto t3 = std::chrono::steady_clock::now();
    double acc = 0.0;
    for (int i = 0; i < bench; ++i)
      acc += sampler->ComputeDepth(queries[i]);
    auto t4 = std::chrono::steady_clock::now();

    // Wave height lookup on the ocean tile.
    auto t5 = std::chrono::steady_clock::now();
    double hacc = 0.0;
    for (int i = 0; i < bench; ++i)
    {
      double h = 0.0;
      wavefield->Height(Eigen::Vector3d(queries[i].x() * 3.0,
          queries[i].y() * 3.0, 0.0), h);
      hacc += h;
    }
    auto t6 = std::chrono::steady_clock::now();

    auto us = [](auto a, auto b) {
      return std::chrono::duration<double, std::micro>(b - a).count();
    };
    out << ",\n  \"bench\": {\n";
    out << "    \"n\": " << bench << ",\n";
    out << "    \"mesh_faces\": " << geom::FaceCount(*sphere) << ",\n";
    out << "    \"ray_mesh_build_us\": " << Num(us(t0, t1)) << ",\n";
    out << "    \"ray_mesh_query_us_per_ray\": " << Num(us(t1, t2) / bench)
        << ",\n";
    out << "    \"ray_mesh_hits\": " << hits << ",\n";
    out << "    \"compute_depth_us_per_query\": " << Num(us(t3, t4) / bench)
        << ",\n";
    out << "    \"wave_height_us_per_query\": " << Num(us(t5, t6) / bench)
        << ",\n";
    out << "    \"checksum\": " << Num(acc + hacc) << "\n";
    out << "  }";
  }

  out << "\n}\n";
  std::cout << out.str();
  return 0;
}
