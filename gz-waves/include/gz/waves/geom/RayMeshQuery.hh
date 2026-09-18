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

/// \file geom/RayMeshQuery.hh
/// \brief Accelerated first-hit ray / mesh intersection.

#ifndef GZ_WAVES_GEOM_RAYMESHQUERY_HH_
#define GZ_WAVES_GEOM_RAYMESHQUERY_HH_

#include <memory>

#include "gz/waves/geom/Types.hh"

namespace gz
{
namespace waves
{
namespace geom
{
class RayMeshQueryPrivate;

/// \brief Acceleration structure over the faces of a mesh answering
/// "first intersection of a ray with the mesh".
///
/// The mesh must outlive the query object.
class RayMeshQuery
{
 public:
  ~RayMeshQuery();

  /// \brief Build the acceleration structure over all faces of _mesh.
  explicit RayMeshQuery(const Mesh& _mesh);

  /// \brief First intersection of the ray (_origin, _direction) with the
  /// mesh. If there is none, the opposite ray is tried, so a point inside
  /// a closed mesh always finds a hit.
  ///
  /// \return true and set _intersection if a hit was found.
  bool FirstIntersection(
      const Point3& _origin,
      const Direction3& _direction,
      Point3& _intersection) const;

 private:
  std::unique_ptr<RayMeshQueryPrivate> data;
};

}  // namespace geom
}  // namespace waves
}  // namespace gz

#endif  // GZ_WAVES_GEOM_RAYMESHQUERY_HH_
