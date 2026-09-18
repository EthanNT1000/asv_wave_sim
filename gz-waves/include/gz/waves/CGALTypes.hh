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

/// \file CGALTypes.hh
/// \brief Compatibility aliases for the historical gz::cgal type names.
///
/// CGAL is no longer used. New code uses gz::waves::geom (geom/Types.hh);
/// these aliases keep the existing public signatures valid.

#ifndef GZ_WAVES_CGALTYPES_HH_
#define GZ_WAVES_CGALTYPES_HH_

#include <memory>

#include "gz/waves/geom/Types.hh"

namespace gz
{
namespace waves
{
namespace geom
{
class RayMeshQuery;
}  // namespace geom
}  // namespace waves

namespace cgal
{
typedef gz::waves::geom::Direction3          Direction3;
typedef gz::waves::geom::Point3              Point3;
typedef gz::waves::geom::Line                Line;
typedef gz::waves::geom::Ray                 Ray;
typedef gz::waves::geom::Triangle            Triangle;
typedef gz::waves::geom::Vector2             Vector2;
typedef gz::waves::geom::Vector3             Vector3;

typedef gz::waves::geom::Mesh                Mesh;
typedef gz::waves::geom::FaceIndex           FaceIndex;
typedef gz::waves::geom::VertexIndex         VertexIndex;
typedef gz::waves::geom::RayMeshQuery        AABBTree;
typedef gz::waves::geom::MeshPtr             MeshPtr;
}  // namespace cgal
}  // namespace gz

#endif  // GZ_WAVES_CGALTYPES_HH_
