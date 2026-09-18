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

#include "gz/waves/MeshTools.hh"

#include <array>
#include <iostream>
#include <iterator>
#include <string>
#include <vector>

#include <gz/common/Console.hh>

#include "gz/waves/geom/Geom.hh"

namespace gz
{
namespace waves
{

//////////////////////////////////////////////////
/// Vertex and Index conventions used by gz::common::Mesh
///
/// Mesh::GetVertexCount()
///   returns the number of vertices in the mesh, which means that
///   the size of _vertices will be 3 * Mesh::GetVertexCount()
///
/// Mesh::GetIndexCount()
///   returns the number of indices in the mesh, which means that
///   the size of _indices will be Mesh::GetIndexCount()
///
void MeshTools::FillArrays(
  const gz::common::Mesh& _source,
  std::vector<float>& _vertices,
  std::vector<int>& _indices
)
{
  double *vertices = nullptr;
  int   *indices  = nullptr;

  // No leaks...
  try
  {
    size_t nv = _source.VertexCount();
    size_t ni = _source.IndexCount();

    // @DEBUG_INFO
    // gzmsg << "nv: " << nv << std::endl;
    // gzmsg << "ni: " << ni << std::endl;

    _source.FillArrays(&vertices, &indices);

    _vertices.insert(_vertices.end(), vertices, vertices + 3 * nv);
    _indices.insert(_indices.end(), indices, indices + ni);
  }
  catch(...)
  {
    gzerr << "Unknown Error in Mesh::FillArrays" << std::endl;
  }
  // Clean up
  if (vertices)
    delete[] vertices;
  if (indices)
    delete[] indices;
}

//////////////////////////////////////////////////
void MeshTools::MakeSurfaceMesh(const gz::common::Mesh& _source,
    geom::Mesh& _target)
{
  std::vector<float> vertices;
  std::vector<int>   indices;

  FillArrays(_source, vertices, indices);

  // Vertices
  for (size_t i=0; i < vertices.size(); )
  {
    auto&& v0 = vertices[i++];
    auto&& v1 = vertices[i++];
    auto&& v2 = vertices[i++];
    geom::Point3 p(v0, v1, v2);

    geom::AddVertex(_target, p);
  }

  // Faces
  for (size_t i=0; i < indices.size(); )
  {
    // @DEBUG_INFO
    // gzmsg << "face" << i/3 << std::endl;

    auto i0 = indices[i++];
    auto i1 = indices[i++];
    auto i2 = indices[i++];
    geom::AddFace(_target, i0, i1, i2);
  }
}

}  // namespace waves
}  // namespace gz
