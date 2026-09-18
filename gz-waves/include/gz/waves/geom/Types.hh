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
/// Phase 2 of the CGAL removal (docs/cgal_audit.md): points and vectors are
/// Eigen (MPL-2.0) fixed-size vectors; Direction3, Line and Triangle are
/// small value types defined here. The mesh container is still a CGAL
/// Surface_mesh until Phase 4; geom/Mesh.hh converts at that boundary.

#ifndef GZ_WAVES_GEOM_TYPES_HH_
#define GZ_WAVES_GEOM_TYPES_HH_

#include <cstdint>
#include <memory>
#include <ostream>

#include <Eigen/Core>

#include <CGAL/Simple_cartesian.h>
#include <CGAL/Surface_mesh.h>

namespace gz
{
namespace waves
{
namespace geom
{
/// \brief Index type for mesh vertices and faces (0-based, contiguous).
typedef int64_t Index;

/// \brief 3D point. Points and vectors share one type, as in most
/// linear-algebra libraries; the names document intent.
typedef Eigen::Vector3d Point3;

/// \brief 3D vector.
typedef Eigen::Vector3d Vector3;

/// \brief 2D vector.
typedef Eigen::Vector2d Vector2;

/// \brief A direction in 3D (a vector whose length carries no meaning).
class Direction3
{
 public:
  Direction3() : v_(0.0, 0.0, 1.0) {}
  Direction3(double x, double y, double z) : v_(x, y, z) {}
  explicit Direction3(const Vector3& v) : v_(v) {}

  /// \brief A vector along the direction (not normalised).
  const Vector3& vector() const { return v_; }

 private:
  Vector3 v_;
};

/// \brief A line through a point with a direction vector.
class Line
{
 public:
  Line() : p_(Vector3::Zero()), v_(0.0, 0.0, 1.0) {}

  /// \brief The line through p and q; point(0) == p, point(1) == q.
  Line(const Point3& p, const Point3& q) : p_(p), v_(q - p) {}

  /// \brief The point at parameter i along the line: p + i * (q - p).
  Point3 point(double i = 0.0) const { return p_ + v_ * i; }

  /// \brief The direction vector q - p.
  const Vector3& to_vector() const { return v_; }

 private:
  Point3 p_;
  Vector3 v_;
};

/// \brief A ray from an origin along a direction.
class Ray
{
 public:
  Ray(const Point3& origin, const Direction3& direction) :
    o_(origin), d_(direction) {}
  const Point3& source() const { return o_; }
  const Direction3& direction() const { return d_; }
  Ray opposite() const { return Ray(o_, Direction3(-d_.vector())); }

 private:
  Point3 o_;
  Direction3 d_;
};

/// \brief A triangle given by three points.
class Triangle
{
 public:
  Triangle() : v_{Vector3::Zero(), Vector3::Zero(), Vector3::Zero()} {}
  Triangle(const Point3& p, const Point3& q, const Point3& r) : v_{p, q, r} {}

  const Point3& operator[](int i) const { return v_[i]; }
  Point3& operator[](int i) { return v_[i]; }
  const Point3& vertex(int i) const { return v_[i]; }

 private:
  Point3 v_[3];
};

inline std::ostream& operator<<(std::ostream& os, const Triangle& t)
{
  return os << "[" << t[0].transpose() << "], [" << t[1].transpose()
            << "], [" << t[2].transpose() << "]";
}

namespace detail
{
/// \brief The CGAL kernel used only by the mesh container and the ray query
/// backend until Phase 4 / Phase 3 replace them.
typedef CGAL::Simple_cartesian<double>  CgalKernel;
typedef CgalKernel::Point_3             CgalPoint3;
}  // namespace detail

/// \brief Indexed triangle mesh (vertex array + faces of three vertices).
typedef CGAL::Surface_mesh<detail::CgalPoint3> Mesh;
typedef std::shared_ptr<Mesh>           MeshPtr;
typedef Mesh::Vertex_index              VertexIndex;
typedef Mesh::Face_index                FaceIndex;

}  // namespace geom
}  // namespace waves
}  // namespace gz

#endif  // GZ_WAVES_GEOM_TYPES_HH_
