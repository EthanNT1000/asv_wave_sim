// Copyright (C) 2019-2026  Rhys Mainwaring and contributors
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

/// \file geom/Mesh.hh
/// \brief Index-based access to a triangle mesh.
///
/// gz-waves only ever treats a mesh as an indexed triangle soup: vertex
/// i has a point, face f has three vertex indices, and indices are stable
/// in 0..N-1 (no deletions). These functions are the only mesh access the
/// library uses. Until Phase 4 the container is a CGAL Surface_mesh storing
/// CGAL points; conversion to the Eigen point type happens here.

#ifndef GZ_WAVES_GEOM_MESH_HH_
#define GZ_WAVES_GEOM_MESH_HH_

#include <array>

#include "gz/waves/geom/Types.hh"

namespace gz
{
namespace waves
{
namespace geom
{
inline Index VertexCount(const Mesh& m)
{
  return static_cast<Index>(m.number_of_vertices());
}

inline Index FaceCount(const Mesh& m)
{
  return static_cast<Index>(m.number_of_faces());
}

inline VertexIndex ToVertexIndex(Index i)
{
  return VertexIndex(static_cast<Mesh::size_type>(i));
}

inline FaceIndex ToFaceIndex(Index i)
{
  return FaceIndex(static_cast<Mesh::size_type>(i));
}

inline Index ToIndex(VertexIndex v)
{
  return static_cast<Index>(static_cast<Mesh::size_type>(v));
}

inline Index ToIndex(FaceIndex f)
{
  return static_cast<Index>(static_cast<Mesh::size_type>(f));
}

inline Point3 VertexPoint(const Mesh& m, Index i)
{
  const detail::CgalPoint3& p = m.point(ToVertexIndex(i));
  return Point3(p.x(), p.y(), p.z());
}

inline void SetVertexPoint(Mesh& m, Index i, const Point3& p)
{
  m.point(ToVertexIndex(i)) = detail::CgalPoint3(p.x(), p.y(), p.z());
}

/// \brief Append a vertex; returns its index.
inline Index AddVertex(Mesh& m, const Point3& p)
{
  return ToIndex(m.add_vertex(detail::CgalPoint3(p.x(), p.y(), p.z())));
}

/// \brief Append a triangular face; returns its index.
inline Index AddFace(Mesh& m, Index i0, Index i1, Index i2)
{
  return ToIndex(m.add_face(
      ToVertexIndex(i0), ToVertexIndex(i1), ToVertexIndex(i2)));
}

/// \brief The three vertex indices of a face, in winding order.
inline std::array<Index, 3> FaceVertices(const Mesh& m, Index f)
{
  std::array<Index, 3> v;
  Mesh::Halfedge_index h = m.halfedge(ToFaceIndex(f));
  v[0] = ToIndex(m.target(h));
  h = m.next(h);
  v[1] = ToIndex(m.target(h));
  h = m.next(h);
  v[2] = ToIndex(m.target(h));
  return v;
}

/// \brief The triangle of a face.
inline Triangle FaceTriangle(const Mesh& m, Index f)
{
  const auto v = FaceVertices(m, f);
  return Triangle(VertexPoint(m, v[0]), VertexPoint(m, v[1]),
      VertexPoint(m, v[2]));
}

}  // namespace geom
}  // namespace waves
}  // namespace gz

#endif  // GZ_WAVES_GEOM_MESH_HH_
