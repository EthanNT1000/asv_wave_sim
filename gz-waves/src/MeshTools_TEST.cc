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

#include <gtest/gtest.h>

#include <chrono>

#include <iostream>
#include <memory>
#include <string>
#include <vector>

#include <gz/common/Console.hh>
#include <gz/common/Util.hh>
#include <gz/common/MeshManager.hh>
#include <gz/common/Mesh.hh>
#include <gz/common/SubMesh.hh>

#include "gz/waves/MeshTools.hh"
#include "gz/waves/Convert.hh"
#include "gz/waves/Geometry.hh"
#include "gz/waves/Grid.hh"
#include "gz/waves/Wavefield.hh"
#include "gz/waves/WaveParameters.hh"
#include "gz/waves/geom/Geom.hh"

namespace geom = gz::waves::geom;

/// \brief Minimal wall-clock timer (seconds).
class Timer
{
 public:
  void start() { start_ = std::chrono::steady_clock::now(); }
  void stop() { stop_ = std::chrono::steady_clock::now(); }
  double time() const
  {
    return std::chrono::duration<double>(stop_ - start_).count();
  }
 private:
  std::chrono::steady_clock::time_point start_, stop_;
};

using gz::waves::Geometry;
using gz::waves::Grid;
using gz::waves::MeshTools;
using gz::waves::ToGz;

//////////////////////////////////////////////////
TEST(MeshTools, FillArraysUnitBox)
{
  // Mesh: 1 x 1 x 1 box
  std::string meshName("box_1x1x1");
  gz::common::MeshManager::Instance()->CreateBox(
    meshName,
    gz::math::Vector3d(1, 1, 1),
    gz::math::Vector2d(1, 1));

  std::vector<float> vertices;
  std::vector<int> indices;
  MeshTools::FillArrays(
    *gz::common::MeshManager::Instance()->MeshByName(meshName),
    vertices, indices);

  // std::cout << "Vertices..." << std::endl;
  // for (auto&& v : vertices)
  //   std::cout << v << std::endl;

  // std::cout << "Indices..." << std::endl;
  // for (auto&& i : indices)
  //   std::cout << i << std::endl;

  // Vertices: 6 sides, 4 points per side, 3 coordinates per point
  EXPECT_EQ(vertices.size(), 6u * 4u * 3u);

  // Indices: 6 sides, 2 faces per side, 3 vertices per face
  EXPECT_EQ(indices.size(), 6u * 2u * 3u);
}

TEST(MeshTools, MakeSurfaceMeshUnitBox)
{
  // Mesh: 1 x 1 x 1 box
  std::string meshName("box_1x1x1");
  gz::common::MeshManager::Instance()->CreateBox(
    meshName,
    gz::math::Vector3d(1, 1, 1),
    gz::math::Vector2d(1, 1));

  geom::Mesh mesh;
  MeshTools::MakeSurfaceMesh(
    *gz::common::MeshManager::Instance()->MeshByName(meshName),
    mesh);

  // std::cout << "Vertices " << std::endl;
  // for(auto&& vertex : mesh.vertices())
  // {
  //   std::cout << vertex << ": " << mesh.point(vertex) << std::endl;
  // }

  // std::cout << "Faces " << std::endl;
  // for(auto&& face : mesh.faces())
  // {
  //   geom::Triangle tri = Geometry::MakeTriangle(mesh, face);
  //   std::cout << face << ": " << tri << std::endl;
  // }

  // 6 sides, 4 vertices per side; 6 sides, 2 triangles per side.
  EXPECT_EQ(geom::VertexCount(mesh), 6 * 4);
  EXPECT_EQ(geom::FaceCount(mesh), 6 * 2);
}

//////////////////////////////////////////////////
void TestExportGridMesh()
{
  std::cout << "TestExportGridMesh..." << std::endl;

  // Wavefield
  Grid grid({500, 500}, {100, 100});

  // Get Mesh
  const auto& mesh = *grid.GetMesh();

  // Create GzMesh
  Timer t;
  t.start();

  std::string name("grid_500m_x_500m_5m_x_5m");
  std::shared_ptr<gz::common::Mesh> gzMesh(new gz::common::Mesh());
  gzMesh->SetName(name);

  std::unique_ptr<gz::common::SubMesh> gzSubMesh(new gz::common::SubMesh());
  int64_t iv = 0;
  const geom::Index nFaces = geom::FaceCount(mesh);
  for (geom::Index face = 0; face < nFaces; ++face)
  {
    geom::Triangle tri  = geom::FaceTriangle(mesh, face);
    geom::Vector3 normal = Geometry::Normal(tri);

    gz::math::Vector3d gzP0(ToGz(tri[0]));
    gz::math::Vector3d gzP1(ToGz(tri[1]));
    gz::math::Vector3d gzP2(ToGz(tri[2]));
    gz::math::Vector3d gzNormal(ToGz(normal));

    gzSubMesh->AddVertex(gzP0);
    gzSubMesh->AddVertex(gzP1);
    gzSubMesh->AddVertex(gzP2);
    gzSubMesh->AddNormal(gzNormal);
    gzSubMesh->AddNormal(gzNormal);
    gzSubMesh->AddNormal(gzNormal);

    gzSubMesh->AddIndex(iv++);
    gzSubMesh->AddIndex(iv++);
    gzSubMesh->AddIndex(iv++);

    // @TODO - calculate texture coordinates
  }

  gzMesh->AddSubMesh(*gzSubMesh);

  t.stop();
  std::cout << "MakeGzMesh: " << t.time() << " sec" << std::endl;

  auto& meshManager = *gz::common::MeshManager::Instance();

  meshManager.Export(
    gzMesh.get(),
    std::string("/Users/rhys/Code/ros/asv_ws/tmp/").append(name),
    "dae");
}

