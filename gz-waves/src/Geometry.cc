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

#include "gz/waves/Geometry.hh"

#include <array>
#include <cmath>
#include <functional>
#include <iostream>
#include <limits>
#include <memory>
#include <string>
#include <variant>

namespace gz
{
namespace waves
{
// Typedefs
double Geometry::TriangleArea(
  const geom::Point3& _p0,
  const geom::Point3& _p1,
  const geom::Point3& _p2
)
{
  geom::Vector3 e1 = _p1 - _p0;
  geom::Vector3 e2 = _p2 - _p0;
  return 0.5 * std::sqrt(geom::SquaredLength(geom::Cross(e1, e2)));
}

double Geometry::TriangleArea(
  const geom::Triangle& _tri
)
{
  return TriangleArea(_tri[0], _tri[1], _tri[2]);
}

geom::Point3 Geometry::TriangleCentroid(
  const geom::Point3& _p0,
  const geom::Point3& _p1,
  const geom::Point3& _p2
)
{
  return geom::Centroid(_p0, _p1, _p2);
}

geom::Point3 Geometry::TriangleCentroid(
  const geom::Triangle& _tri
)
{
  return geom::Centroid(_tri[0], _tri[1], _tri[2]);
}

geom::Point3 Geometry::MidPoint(
  const geom::Point3& _p0,
  const geom::Point3& _p1
)
{
  return geom::MidPoint(_p0, _p1);
}

geom::Vector2 Geometry::Normalize(const geom::Vector2& _v)
{
  if (_v == geom::NullVector2())
    return _v;
  else
    return _v/std::sqrt(geom::SquaredLength(_v));
}

geom::Vector3 Geometry::Normalize(const geom::Vector3& _v)
{
  if (_v == geom::NullVector())
    return _v;
  else
    return _v/std::sqrt(geom::SquaredLength(_v));
}

namespace
{
/// \brief Unit normal of the triangle (p0, p1, p2), or the null vector if
/// the points are degenerate.
///
/// A triangle is degenerate when geom::Collinear reports it: the cross
/// product of its edge vectors is exactly zero (the CGAL Cartesian
/// predicate on doubles). Near-degenerate triangles, whose edge vectors
/// are almost but not exactly parallel, keep the normal of that cross
/// product, as they did with CGAL. No tolerance is applied, so clipped
/// sub-triangles are classified exactly as before.
/// See https://github.com/srmainwaring/asv_wave_sim/issues/50.
geom::Vector3 UnitNormalOrNull(
  const geom::Point3& _p0,
  const geom::Point3& _p1,
  const geom::Point3& _p2)
{
  if (geom::Collinear(_p0, _p1, _p2))
    return geom::NullVector();
  auto n = geom::UnitlessNormal(_p0, _p1, _p2);
  return n/std::sqrt(geom::SquaredLength(n));
}
}  // namespace

geom::Vector3 Geometry::Normal(
  const geom::Point3& _p0,
  const geom::Point3& _p1,
  const geom::Point3& _p2
)
{
  return UnitNormalOrNull(_p0, _p1, _p2);
}

geom::Vector3 Geometry::Normal(
  const geom::Triangle& _tri
)
{
  return UnitNormalOrNull(_tri[0], _tri[1], _tri[2]);
}

geom::Vector3 Geometry::Normal(const geom::Mesh& _mesh, geom::FaceIndex _face)
{
  const auto v = geom::FaceVertices(_mesh, geom::ToIndex(_face));
  return UnitNormalOrNull(
      geom::VertexPoint(_mesh, v[0]),
      geom::VertexPoint(_mesh, v[1]),
      geom::VertexPoint(_mesh, v[2]));
}

geom::Point3 Geometry::HorizontalIntercept(
  const geom::Point3& _high,
  const geom::Point3& _mid,
  const geom::Point3& _low
)
{
  // If _high.Z() = _low().Z() then _high.Z() = _low.Z() = _mid.Z()
  // so any point on the line LH will satisfy the intercept condition,
  // including t=0, which we set as default (and so return _low).
  double t = 0.0;
  double div = _high.z() - _low.z();
  if (std::abs(div) > std::numeric_limits<double>::epsilon())
  {
    t = (_mid.z() - _low.z()) / div;
  }
  return geom::Point3(
    _low.x() + t * (_high.x() - _low.x()),
    _low.y() + t * (_high.y() - _low.y()),
    _mid.z());
}

bool Geometry::RayIntersectsTriangle(
  const geom::Point3& _origin,
  const geom::Direction3& _direction,
  const geom::Triangle& _tri,
  geom::Point3& _intersection
)
{
  const geom::Point3& p0 = _tri[0];
  const geom::Point3& p1 = _tri[1];
  const geom::Point3& p2 = _tri[2];
  geom::Vector3 _ray(_direction.vector());
  geom::Vector3 e1 = p1 - p0;
  geom::Vector3 e2 = p2 - p0;
  geom::Vector3 h = geom::Cross(_ray, e2);
  double a = geom::Dot(e1, h);
  if (a > -std::numeric_limits<double>::epsilon()
    && a < std::numeric_limits<double>::epsilon())
      return false;    // This ray is parallel to this triangle.
  double f = 1.0 / a;
  geom::Vector3 s = _origin - p0;
  double u = f * geom::Dot(s, h);
  if (u < 0.0 || u > 1.0)
      return false;
  geom::Vector3 q = geom::Cross(s, e1);
  double v = f * geom::Dot(_ray, q);
  if (v < 0.0 || u + v > 1.0)
      return false;
  // At this stage we can compute t to find out where the intersection point
  // is on the line.
  double t = f * geom::Dot(e2, q);
  if (t > std::numeric_limits<double>::epsilon())
  {
    // ray intersection
    _intersection = _origin + _ray * t;
    return true;
  } else {
    // This means that there is a line intersection but not a ray intersection.
    return false;
  }
}

bool Geometry::LineIntersectsTriangle(
  const geom::Point3& _origin,
  const geom::Direction3& _direction,
  const geom::Triangle& _tri,
  geom::Point3& _intersection
)
{
  return LineIntersectsTriangle(
    _origin, _direction, _tri[0], _tri[1], _tri[2], _intersection);
}

bool Geometry::LineIntersectsTriangle(
  const geom::Point3& _origin,
  const geom::Direction3& _direction,
  const geom::Point3& _p0,
  const geom::Point3& _p1,
  const geom::Point3& _p2,
  geom::Point3& _intersection
)
{
  geom::Vector3 _ray(_direction.vector());
  geom::Vector3 e1 = _p1 - _p0;
  geom::Vector3 e2 = _p2 - _p0;
  geom::Vector3 h = geom::Cross(_ray, e2);
  double a = geom::Dot(e1, h);
  if (a > -std::numeric_limits<double>::epsilon()
    && a < std::numeric_limits<double>::epsilon())
      return false;    // This ray is parallel to this triangle.
  double f = 1.0 / a;
  geom::Vector3 s = _origin - _p0;
  double u = f * geom::Dot(s, h);
  if (u < 0.0 || u > 1.0)
      return false;
  geom::Vector3 q = geom::Cross(s, e1);
  double v = f * geom::Dot(_ray, q);
  if (v < 0.0 || u + v > 1.0)
      return false;
  // At this stage we can compute t to find out where the intersection point
  // is on the line.
  double t = f * geom::Dot(e2, q);
  _intersection = _origin + _ray * t;
  return true;
}

geom::Triangle Geometry::MakeTriangle(
    const geom::Mesh& _mesh, geom::FaceIndex _face)
{
  return geom::FaceTriangle(_mesh, geom::ToIndex(_face));
}

std::shared_ptr<geom::RayMeshQuery> Geometry::MakeAABBTree(
    const geom::Mesh& _mesh)
{
  return std::make_shared<geom::RayMeshQuery>(_mesh);
}

//////////////////////////////////////////////////
bool Geometry::SearchMesh(
  const geom::Mesh& _mesh,
  const geom::Point3& _origin,
  const geom::Direction3& _direction,
  geom::Point3& _intersection
)
{
  geom::RayMeshQuery query(_mesh);
  return query.FirstIntersection(_origin, _direction, _intersection);
}

//////////////////////////////////////////////////
bool Geometry::SearchMesh(
  const geom::RayMeshQuery& _tree,
  const geom::Point3& _origin,
  const geom::Direction3& _direction,
  geom::Point3& _intersection
)
{
  return _tree.FirstIntersection(_origin, _direction, _intersection);
}

}  // namespace waves
}  // namespace gz
