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

/// \file geom/Types.hh
/// \brief Geometry kernel types used throughout gz-waves.
///
/// Every geometric type the library uses is declared here under
/// gz::waves::geom. Code outside geom:: must not name the backing
/// library directly; it uses these aliases and the free functions in
/// geom/Vector.hh, geom/Mesh.hh and geom/RayMeshQuery.hh.
///
/// Phase 0 of the CGAL removal (docs/cgal_audit.md): the aliases resolve
/// to the CGAL Simple_cartesian<double> kernel so behaviour is unchanged.

#ifndef GZ_WAVES_GEOM_TYPES_HH_
#define GZ_WAVES_GEOM_TYPES_HH_

#include <cstdint>

#include "gz/waves/CGALTypes.hh"

namespace gz
{
namespace waves
{
namespace geom
{
/// \brief Index type for mesh vertices and faces (0-based, contiguous).
typedef int64_t Index;

typedef cgal::Point3      Point3;
typedef cgal::Vector3     Vector3;
typedef cgal::Vector2     Vector2;
typedef cgal::Direction3  Direction3;
typedef cgal::Line        Line;
typedef cgal::Ray         Ray;
typedef cgal::Triangle    Triangle;

/// \brief Indexed triangle mesh (vertex array + faces of three vertices).
typedef cgal::Mesh        Mesh;
typedef cgal::MeshPtr     MeshPtr;
typedef cgal::VertexIndex VertexIndex;
typedef cgal::FaceIndex   FaceIndex;

}  // namespace geom
}  // namespace waves
}  // namespace gz

#endif  // GZ_WAVES_GEOM_TYPES_HH_
