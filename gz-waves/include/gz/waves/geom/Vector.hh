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

/// \file geom/Vector.hh
/// \brief Vector algebra and elementary triangle functions on geom types.
///
/// The formulas match the CGAL Simple_cartesian<double> operations they
/// replaced term by term (cross product component order, (p+q+r)/3
/// centroid, (p+q)/2 midpoint, exact-zero collinearity test), so results are
/// unchanged to the last bit where the compiler evaluates them in the same
/// order.

#ifndef GZ_WAVES_GEOM_VECTOR_HH_
#define GZ_WAVES_GEOM_VECTOR_HH_

#include <cmath>

#include "gz/waves/geom/Types.hh"

namespace gz
{
namespace waves
{
namespace geom
{
/// \brief The origin (0, 0, 0) as a point.
inline Point3 Origin() { return Point3::Zero(); }

/// \brief The zero vector.
inline Vector3 NullVector() { return Vector3::Zero(); }

/// \brief The zero 2D vector.
inline Vector2 NullVector2() { return Vector2::Zero(); }

/// \brief Convert a kernel scalar to double (identity for double kernels).
inline double ToDouble(double x) { return x; }

inline Vector3 Cross(const Vector3& a, const Vector3& b)
{
  return Vector3(a.y() * b.z() - a.z() * b.y(),
                 a.z() * b.x() - a.x() * b.z(),
                 a.x() * b.y() - a.y() * b.x());
}

inline double Dot(const Vector3& a, const Vector3& b)
{
  return a.x() * b.x() + a.y() * b.y() + a.z() * b.z();
}

inline double SquaredLength(const Vector3& v)
{
  return Dot(v, v);
}

inline double SquaredLength(const Vector2& v)
{
  return v.x() * v.x() + v.y() * v.y();
}

inline double Length(const Vector3& v)
{
  return std::sqrt(SquaredLength(v));
}

/// \brief Point at the given coordinates.
inline Point3 MakePoint(double x, double y, double z)
{
  return Point3(x, y, z);
}

inline Vector3 MakeVector(double x, double y, double z)
{
  return Vector3(x, y, z);
}

/// \brief Direction from a vector (need not be unit length).
inline Direction3 MakeDirection(double x, double y, double z)
{
  return Direction3(x, y, z);
}

/// \brief The vector representation of a direction.
inline const Vector3& DirectionVector(const Direction3& d)
{
  return d.vector();
}

/// \brief Unnormalised normal (q - p) x (r - p).
inline Vector3 UnitlessNormal(const Point3& p, const Point3& q,
    const Point3& r)
{
  return Cross(q - p, r - p);
}

/// \brief True if the three points are exactly collinear: all three 2x2
/// minors of (q - p, r - p) vanish (the CGAL Cartesian predicate).
inline bool Collinear(const Point3& p, const Point3& q, const Point3& r)
{
  const Vector3 n = Cross(q - p, r - p);
  return n.x() == 0.0 && n.y() == 0.0 && n.z() == 0.0;
}

inline Point3 Centroid(const Point3& p, const Point3& q, const Point3& r)
{
  return Point3((p.x() + q.x() + r.x()) / 3.0,
                (p.y() + q.y() + r.y()) / 3.0,
                (p.z() + q.z() + r.z()) / 3.0);
}

inline Point3 MidPoint(const Point3& p, const Point3& q)
{
  return Point3((p.x() + q.x()) / 2.0,
                (p.y() + q.y()) / 2.0,
                (p.z() + q.z()) / 2.0);
}

/// \brief Construct a line through two points.
inline Line MakeLine(const Point3& p, const Point3& q)
{
  return Line(p, q);
}

/// \brief A point on the line.
inline Point3 LinePoint(const Line& l) { return l.point(); }

/// \brief The direction vector of the line (q - p for MakeLine(p, q)).
inline const Vector3& LineVector(const Line& l) { return l.to_vector(); }

}  // namespace geom
}  // namespace waves
}  // namespace gz

#endif  // GZ_WAVES_GEOM_VECTOR_HH_
