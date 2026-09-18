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
/// library uses.

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
inline Index VertexCount(const Mesh& m) { return m.VertexCount(); }

inline Index FaceCount(const Mesh& m) { return m.FaceCount(); }

inline VertexIndex ToVertexIndex(Index i) { return i; }

inline FaceIndex ToFaceIndex(Index i) { return i; }

inline Index ToIndex(Index i) { return i; }

inline const Point3& VertexPoint(const Mesh& m, Index i)
{
  return m.Point(i);
}

inline void SetVertexPoint(Mesh& m, Index i, const Point3& p)
{
  m.Point(i) = p;
}

/// \brief Append a vertex; returns its index.
inline Index AddVertex(Mesh& m, const Point3& p) { return m.AddVertex(p); }

/// \brief Append a triangular face; returns its index.
inline Index AddFace(Mesh& m, Index i0, Index i1, Index i2)
{
  return m.AddFace(i0, i1, i2);
}

/// \brief The three vertex indices of a face, in winding order.
inline const std::array<Index, 3>& FaceVertices(const Mesh& m, Index f)
{
  return m.FaceAt(f);
}

/// \brief The triangle of a face.
inline Triangle FaceTriangle(const Mesh& m, Index f)
{
  const auto& v = m.FaceAt(f);
  return Triangle(m.Point(v[0]), m.Point(v[1]), m.Point(v[2]));
}

}  // namespace geom
}  // namespace waves
}  // namespace gz

#endif  // GZ_WAVES_GEOM_MESH_HH_
