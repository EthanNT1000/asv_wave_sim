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

/// \file CGALTypes.hh
/// \brief Type definitions for CGAL structures used in the library.

#ifndef GZ_WAVES_CGALTYPES_HH_
#define GZ_WAVES_CGALTYPES_HH_

#include <CGAL/Simple_cartesian.h>
#include <CGAL/Surface_mesh.h>

#include <memory>

namespace gz
{
namespace waves
{
namespace geom
{
class RayMeshQuery;
}  // namespace geom
}  // namespace waves

/// \brief Compatibility aliases. New code uses gz::waves::geom (see
/// geom/Types.hh); these names are kept so existing signatures compile.
namespace cgal
{
//////////////////////////////////////////////////
// CGAL Typedefs

// 2D/3D Linear Geometry
typedef CGAL::Simple_cartesian<double>  Kernel;
typedef Kernel::Direction_2             Direction2;
typedef Kernel::Direction_3             Direction3;
typedef Kernel::Point_3                 Point3;
typedef Kernel::Line_3                  Line;
typedef Kernel::Ray_3                   Ray;
typedef Kernel::Triangle_3              Triangle;
typedef Kernel::Vector_2                Vector2;
typedef Kernel::Vector_3                Vector3;

// SurfaceMesh
typedef CGAL::Surface_mesh<Point3>      Mesh;
typedef Mesh::Face_index                FaceIndex;
typedef Mesh::Halfedge_index            HalfedgeIndex;
typedef Mesh::Vertex_index              VertexIndex;

// AABB Tree
typedef gz::waves::geom::RayMeshQuery AABBTree;

// Pointers
typedef std::shared_ptr<Mesh>           MeshPtr;
}  // namespace cgal
}  // namespace gz

#endif  // GZ_WAVES_CGALTYPES_HH_
