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

/// \file geom/RayMeshQuery.cc
/// \brief Embree backend for first-hit ray / mesh queries (Phase 3 of the
/// CGAL removal, docs/cgal_audit.md).
///
/// Embree (Apache-2.0) builds the bounding-volume hierarchy and finds the
/// nearest hit primitive in single precision. The intersection point is then
/// recomputed in double precision on that triangle (ray / plane
/// intersection), so the result matches the previous double-precision CGAL
/// AABB-tree answer to rounding. Vertices are stored relative to the mesh
/// bounding-box centre to keep the float traversal accurate for meshes far
/// from the origin.

#include "gz/waves/geom/RayMeshQuery.hh"

#include <cmath>
#include <limits>
#include <vector>

#if __has_include(<embree4/rtcore.h>)
#include <embree4/rtcore.h>
#else
#include <embree3/rtcore.h>
#endif

#include "gz/waves/geom/Mesh.hh"
#include "gz/waves/geom/Vector.hh"

namespace gz
{
namespace waves
{
namespace geom
{
namespace
{
/// \brief One Embree device shared by all queries in the process.
RTCDevice SharedDevice()
{
  static RTCDevice device = rtcNewDevice(nullptr);
  return device;
}

/// \brief Ray / plane intersection with the plane of triangle (a, b, c),
/// in double precision. Returns false if the ray is parallel to the plane.
bool RefineHit(const Point3& o, const Vector3& d,
    const Point3& a, const Point3& b, const Point3& c, Point3& hit)
{
  const Vector3 n = Cross(b - a, c - a);
  const double denom = Dot(d, n);
  if (denom == 0.0)
    return false;
  const double t = Dot(a - o, n) / denom;
  hit = o + d * t;
  return true;
}
}  // namespace

class RayMeshQueryPrivate
{
 public:
  explicit RayMeshQueryPrivate(const Mesh& _mesh) :
    mesh(_mesh), scene(nullptr), center(Vector3::Zero())
  {
    const Index nV = VertexCount(_mesh);
    const Index nF = FaceCount(_mesh);
    if (nV == 0 || nF == 0)
      return;

    // Bounding-box centre, used as the float origin.
    Point3 lo = VertexPoint(_mesh, 0), hi = lo;
    for (Index i = 1; i < nV; ++i)
    {
      const Point3 p = VertexPoint(_mesh, i);
      lo = lo.cwiseMin(p);
      hi = hi.cwiseMax(p);
    }
    center = 0.5 * (lo + hi);

    RTCDevice device = SharedDevice();
    scene = rtcNewScene(device);
    RTCGeometry geometry = rtcNewGeometry(device, RTC_GEOMETRY_TYPE_TRIANGLE);

    float* vb = static_cast<float*>(rtcSetNewGeometryBuffer(geometry,
        RTC_BUFFER_TYPE_VERTEX, 0, RTC_FORMAT_FLOAT3, 3 * sizeof(float),
        static_cast<size_t>(nV)));
    for (Index i = 0; i < nV; ++i)
    {
      const Point3 p = VertexPoint(_mesh, i) - center;
      vb[3 * i + 0] = static_cast<float>(p.x());
      vb[3 * i + 1] = static_cast<float>(p.y());
      vb[3 * i + 2] = static_cast<float>(p.z());
    }

    unsigned* ib = static_cast<unsigned*>(rtcSetNewGeometryBuffer(geometry,
        RTC_BUFFER_TYPE_INDEX, 0, RTC_FORMAT_UINT3, 3 * sizeof(unsigned),
        static_cast<size_t>(nF)));
    faces.reserve(nF);
    for (Index f = 0; f < nF; ++f)
    {
      const auto& v = FaceVertices(_mesh, f);
      ib[3 * f + 0] = static_cast<unsigned>(v[0]);
      ib[3 * f + 1] = static_cast<unsigned>(v[1]);
      ib[3 * f + 2] = static_cast<unsigned>(v[2]);
      faces.push_back(v);
    }

    rtcCommitGeometry(geometry);
    rtcAttachGeometry(scene, geometry);
    rtcReleaseGeometry(geometry);
    rtcCommitScene(scene);
  }

  ~RayMeshQueryPrivate()
  {
    if (scene)
      rtcReleaseScene(scene);
  }

  /// \brief Nearest hit along the ray (o, d); false if none.
  bool Intersect(const Point3& o, const Vector3& d, Point3& hit) const
  {
    if (!scene)
      return false;

    const Point3 ol = o - center;
    RTCRayHit rh;
    rh.ray.org_x = static_cast<float>(ol.x());
    rh.ray.org_y = static_cast<float>(ol.y());
    rh.ray.org_z = static_cast<float>(ol.z());
    rh.ray.dir_x = static_cast<float>(d.x());
    rh.ray.dir_y = static_cast<float>(d.y());
    rh.ray.dir_z = static_cast<float>(d.z());
    rh.ray.tnear = 0.0f;
    rh.ray.tfar = std::numeric_limits<float>::infinity();
    rh.ray.time = 0.0f;
    rh.ray.mask = static_cast<unsigned>(-1);
    rh.ray.id = 0;
    rh.ray.flags = 0;
    rh.hit.geomID = RTC_INVALID_GEOMETRY_ID;
    rh.hit.primID = RTC_INVALID_GEOMETRY_ID;
    rh.hit.instID[0] = RTC_INVALID_GEOMETRY_ID;

#if RTC_VERSION_MAJOR >= 4
    rtcIntersect1(scene, &rh);
#else
    RTCIntersectContext context;
    rtcInitIntersectContext(&context);
    rtcIntersect1(scene, &context, &rh);
#endif

    if (rh.hit.geomID == RTC_INVALID_GEOMETRY_ID)
      return false;

    // Recompute the hit in double on the reported triangle.
    const auto& v = faces[rh.hit.primID];
    if (RefineHit(o, d, VertexPoint(mesh, v[0]), VertexPoint(mesh, v[1]),
        VertexPoint(mesh, v[2]), hit))
    {
      return true;
    }
    hit = o + d * static_cast<double>(rh.ray.tfar);
    return true;
  }

  const Mesh& mesh;
  RTCScene scene;
  Vector3 center;
  std::vector<std::array<Index, 3>> faces;
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
  const Vector3& d = _direction.vector();
  if (data->Intersect(_origin, d, _intersection))
    return true;
  // Search the opposite direction as well.
  return data->Intersect(_origin, Vector3(-d), _intersection);
}

}  // namespace geom
}  // namespace waves
}  // namespace gz
