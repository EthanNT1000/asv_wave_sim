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

#include "gz/waves/geom/RayMeshQuery.hh"

#include <CGAL/AABB_face_graph_triangle_primitive.h>
#if CGAL_VERSION_MAJOR >= 6
#include <CGAL/AABB_traits_3.h>
#include <optional>
#include <variant>
#else
#include <CGAL/AABB_traits.h>
#include <boost/optional.hpp>
#include <boost/variant.hpp>
#endif
#include <CGAL/AABB_tree.h>

namespace gz
{
namespace waves
{
namespace geom
{
// Phase 0 backend: the CGAL AABB tree (docs/cgal_audit.md, Sec. 1a).
typedef CGAL::AABB_face_graph_triangle_primitive<cgal::Mesh> Primitive;
#if CGAL_VERSION_MAJOR >= 6
typedef CGAL::AABB_traits_3<cgal::Kernel, Primitive> Traits;
#else
typedef CGAL::AABB_traits<cgal::Kernel, Primitive> Traits;
#endif
typedef CGAL::AABB_tree<Traits> Tree;

class RayMeshQueryPrivate
{
 public:
  explicit RayMeshQueryPrivate(const Mesh& _mesh) :
    tree(faces(_mesh).first, faces(_mesh).second, _mesh)
  {
    tree.build();
  }

  Tree tree;
};

//////////////////////////////////////////////////
RayMeshQuery::~RayMeshQuery() = default;

//////////////////////////////////////////////////
RayMeshQuery::RayMeshQuery(const Mesh& _mesh) :
  data(new RayMeshQueryPrivate(_mesh))
{
}

//////////////////////////////////////////////////
bool RayMeshQuery::FirstIntersection(
    const Point3& _origin,
    const Direction3& _direction,
    Point3& _intersection) const
{
#if CGAL_VERSION_MAJOR >= 6
  typedef std::optional<Tree::Intersection_and_primitive_id<
      cgal::Ray>::Type> RayIntersection;
#else
  typedef boost::optional<Tree::Intersection_and_primitive_id<
      cgal::Ray>::Type> RayIntersection;
#endif

  cgal::Ray query(_origin, _direction);
  RayIntersection intersection = data->tree.first_intersection(query);

  // Search both directions
  if (!intersection)
    intersection = data->tree.first_intersection(query.opposite());

  if (intersection)
  {
#if CGAL_VERSION_MAJOR >= 6
    const cgal::Point3* p = std::get_if<cgal::Point3>(&(intersection->first));
#else
    const cgal::Point3* p = boost::get<cgal::Point3>(&(intersection->first));
#endif
    if (p)
    {
      _intersection = *p;
      return true;
    }
  }
  return false;
}

}  // namespace geom
}  // namespace waves
}  // namespace gz
