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

#include "gz/waves/TriangulatedGrid.hh"

#include <algorithm>
#include <cmath>
#include <iostream>
#include <memory>
#include <utility>
#include <vector>

#include "gz/waves/Geometry.hh"
#include "gz/waves/geom/Geom.hh"

namespace gz
{
namespace waves
{

/// \internal
/// Regular-lattice implementation (Phase 1 of the CGAL removal, see
/// docs/cgal_audit.md Sec. 3).
///
/// The tile is an (nx+1) x (ny+1) lattice with cell size (lx/nx, ly/ny) whose
/// cells are split by the idx0-idx2 diagonal. Point location is index
/// arithmetic on the undisplaced lattice followed by a ray/triangle test on
/// the *current* vertex positions, expanding cell by cell in a ring around the
/// guessed cell. The ring search is required because the wave simulations
/// displace the vertices horizontally (Gerstner / FFT choppiness), so the
/// triangle containing a query may be the neighbour of the lattice cell that
/// contains it. The located triangle and the height computed from it are the
/// same as the constrained Delaunay triangulation previously used, because the
/// triangulation faces were exactly the lattice cells split by the constrained
/// diagonal.
class TriangulatedGrid::Private {
 public:
  ~Private();
  Private(Index nx, Index ny, double lx, double ly);
  void CreateMesh();
  void CreateTriangulation();

  bool Locate(const geom::Point3& query, int64_t& faceIndex) const;
  bool Height(const geom::Point3& query, double& height) const;
  bool Height(const std::vector<geom::Point3>& queries,
      std::vector<double>& heights) const;
  bool Interpolate(TriangulatedGrid& patch) const;

  const Point3Range& Points() const;
  const Index3Range& Indices() const;
  const geom::Point3& Origin() const;
  void ApplyPose(const gz::math::Pose3d& pose);

  bool IsValid(bool verbose = false) const;
  void DebugPrintMesh() const;
  void DebugPrintTriangulation() const;
  void UpdatePoints(const std::vector<geom::Point3>& points);
  void UpdatePoints(const std::vector<gz::math::Vector3d>& from);
  void UpdatePoints(const geom::Mesh& from);

  /// \brief Locate the triangle containing the xy of query and intersect
  /// the vertical line through query with it.
  ///
  /// \param[in] query      Query point.
  /// \param[in, out] hint  Face index to try first (-1 for none); set to the
  ///                       found face on success.
  /// \param[out] intersection The point on the surface above/below query.
  /// \return true if a containing triangle was found.
  bool LocateAndIntersect(const geom::Point3& query, int64_t& hint,
      geom::Point3& intersection) const;

  /// \brief Test the two triangles of cell (ix, iy).
  bool IntersectCell(Index ix, Index iy, const geom::Point3& query,
      int64_t& faceIndex, geom::Point3& intersection) const;

  /// \brief Test one triangle.
  bool IntersectFace(int64_t faceIndex, const geom::Point3& query,
      geom::Point3& intersection) const;

  // Dimensions
  Index nx_;
  Index ny_;
  double lx_;
  double ly_;

  // Mesh
  geom::Point3              origin_;
  std::vector<geom::Point3> points0_;
  std::vector<geom::Point3> points_;
  std::vector<Index3> indices_;
  std::vector<Index3> infinite_indices_;
};

TriangulatedGrid::Private::~Private() {
}

TriangulatedGrid::Private::Private(Index nx, Index ny, double lx, double ly) :
  nx_(nx), ny_(ny), lx_(lx), ly_(ly), origin_(geom::Origin()) {
}

void TriangulatedGrid::Private::CreateMesh() {
  double dlx = lx_ / nx_;
  double dly = ly_ / ny_;
  double lxm = - lx_ / 2.0;
  double lym = - ly_ / 2.0;
  const Index nx_plus1 = nx_ + 1;

  // Points - (nx_+1) points_ in each row / column
  for (int64_t iy=0; iy <= ny_; ++iy) {
    double py = iy * dly + lym;
    for (int64_t ix=0; ix <= nx_; ++ix) {
      // Vertex position
      double px = ix * dlx + lxm;
      geom::Point3 point(px, py, 0.0);
      points_.push_back(point);
    }
  }
  // Copy initial points
  points0_ = points_;

  // Face indices
  for (int64_t iy=0; iy < ny_; ++iy) {
    for (int64_t ix=0; ix < nx_; ++ix) {
      // Get the points in the cell coordinates
      int64_t idx0 = iy * nx_plus1 + ix;
      int64_t idx1 = iy * nx_plus1 + ix + 1;
      int64_t idx2 = (iy+1) * nx_plus1 + ix + 1;
      int64_t idx3 = (iy+1) * nx_plus1 + ix;

      // Face indices
      indices_.push_back({ idx0, idx1, idx2 });
      indices_.push_back({ idx0, idx2, idx3 });
    }
  }

  // Infinite indices (follow edges counter clockwise around grid)
  infinite_indices_.resize(2*nx_ + 2*ny_);
  for (int64_t ix=0; ix < nx_; ++ix) {
    // bottom (0..nx-1)
    int64_t idx = ix;
    infinite_indices_[ix] = { idx, idx + 1 };

    // top (nx+ny..2*nx+ny-1)
    idx = ny_ * nx_plus1 + ix;
    infinite_indices_[2*nx_ + ny_ - 1 - ix] = { idx + 1, idx };
  }
  for (int64_t iy=0; iy < ny_; ++iy) {
    // right (nx..nx+ny-1)
    int64_t idx = iy * nx_plus1 + nx_;
    infinite_indices_[nx_ + iy] = { idx, idx + nx_plus1 };

    // left (2*nx+ny..2*nx+2*ny-1)
    idx = iy * nx_plus1;
    infinite_indices_[2*nx_ + 2*ny_ - 1 - iy] = { idx + nx_plus1, idx };
  }
}

void TriangulatedGrid::Private::CreateTriangulation() {
  // The lattice connectivity built in CreateMesh is the triangulation:
  // face 2*(iy*nx + ix) + k is triangle k of cell (ix, iy). Nothing to build.
}

bool TriangulatedGrid::Private::IntersectFace(int64_t faceIndex,
    const geom::Point3& query, geom::Point3& intersection) const {
  const Index3& f = indices_[faceIndex];
  const geom::Direction3 direction(0, 0, 1);
  return Geometry::LineIntersectsTriangle(query, direction,
      points_[f[0]], points_[f[1]], points_[f[2]], intersection);
}

bool TriangulatedGrid::Private::IntersectCell(Index ix, Index iy,
    const geom::Point3& query, int64_t& faceIndex,
    geom::Point3& intersection) const {
  const int64_t f0 = 2 * (iy * nx_ + ix);
  if (IntersectFace(f0, query, intersection)) {
    faceIndex = f0;
    return true;
  }
  if (IntersectFace(f0 + 1, query, intersection)) {
    faceIndex = f0 + 1;
    return true;
  }
  return false;
}

bool TriangulatedGrid::Private::LocateAndIntersect(const geom::Point3& query,
    int64_t& hint, geom::Point3& intersection) const {
  if (indices_.empty())
    return false;

  // Hint from the previous query (spatially coherent queries).
  if (hint >= 0 && hint < static_cast<int64_t>(indices_.size()) &&
      IntersectFace(hint, query, intersection)) {
    return true;
  }

  // Index math on the undisplaced lattice: cell containing the query xy.
  const double dlx = lx_ / nx_;
  const double dly = ly_ / ny_;
  const double x0 = origin_.x() - lx_ / 2.0;
  const double y0 = origin_.y() - ly_ / 2.0;
  Index ix = static_cast<Index>(std::floor((query.x() - x0) / dlx));
  Index iy = static_cast<Index>(std::floor((query.y() - y0) / dly));
  ix = std::min(std::max<Index>(ix, 0), nx_ - 1);
  iy = std::min(std::max<Index>(iy, 0), ny_ - 1);

  int64_t face = -1;
  if (IntersectCell(ix, iy, query, face, intersection)) {
    hint = face;
    return true;
  }

  // Expand in rings around the guessed cell. The horizontal displacement of
  // the wave vertices is bounded by the wave steepness, so in practice the
  // containing triangle is found within one or two rings; the loop
  // nevertheless covers the whole tile so a query is never missed while some
  // triangle contains it.
  const Index maxRing = std::max(nx_, ny_);
  for (Index r = 1; r <= maxRing; ++r) {
    const Index ixm = ix - r, ixp = ix + r, iym = iy - r, iyp = iy + r;
    if (ixm < 0 && ixp >= nx_ && iym < 0 && iyp >= ny_)
      break;
    // bottom and top rows of the ring
    const Index i0 = std::max<Index>(ixm, 0);
    const Index i1 = std::min<Index>(ixp, nx_ - 1);
    for (Index i = i0; i <= i1; ++i) {
      if (iym >= 0 && IntersectCell(i, iym, query, face, intersection)) {
        hint = face; return true;
      }
      if (iyp < ny_ && IntersectCell(i, iyp, query, face, intersection)) {
        hint = face; return true;
      }
    }
    // left and right columns (excluding corners already visited)
    const Index j0 = std::max<Index>(iym + 1, 0);
    const Index j1 = std::min<Index>(iyp - 1, ny_ - 1);
    for (Index j = j0; j <= j1; ++j) {
      if (ixm >= 0 && IntersectCell(ixm, j, query, face, intersection)) {
        hint = face; return true;
      }
      if (ixp < nx_ && IntersectCell(ixp, j, query, face, intersection)) {
        hint = face; return true;
      }
    }
  }
  return false;
}

bool TriangulatedGrid::Private::Locate(
    const geom::Point3& query, int64_t& faceIndex) const {
  int64_t hint = -1;
  geom::Point3 intersection(query);
  if (LocateAndIntersect(query, hint, intersection)) {
    faceIndex = hint;
    return true;
  }
  return false;
}

bool TriangulatedGrid::Private::Height(
    const geom::Point3& query, double& height) const {
  height = 0.0;
  int64_t hint = -1;
  geom::Point3 intersection(query);
  bool found = LocateAndIntersect(query, hint, intersection);
  if (found) {
    height = intersection.z() - query.z();
  }
  return found;
}

bool TriangulatedGrid::Private::Height(
    const std::vector<geom::Point3>& queries,
    std::vector<double>& heights) const
{
  bool foundAll = true;
  int64_t hint = -1;
  for (uint64_t i=0; i < heights.size(); ++i)
  {
    double height_i = 0.0;
    const geom::Point3& query = queries[i];
    geom::Point3 intersection(query);
    bool found = LocateAndIntersect(query, hint, intersection);
    if (found) {
      height_i = intersection.z() - query.z();
    }
    heights[i] = height_i;
    foundAll &= found;
  }
  return foundAll;
}

bool TriangulatedGrid::Private::Interpolate(TriangulatedGrid& patch) const {
  bool foundAll = true;

  int64_t hint = -1;
  for (auto it = patch.impl_->points_.begin();
      it != patch.impl_->points_.end(); ++it) {
    double height = 0.0;
    const geom::Point3& query = *it;
    geom::Point3 intersection(query);
    bool found = LocateAndIntersect(query, hint, intersection);
    if (found) {
      height = intersection.z();
    }
    // @NOTE this assumes the patch initially has height = 0.0;
    *it = geom::Point3(query.x(), query.y(), height);

    foundAll &= found;
  }
  return foundAll;
}

const Point3Range& TriangulatedGrid::Private::Points() const {
  return points_;
}

const Index3Range& TriangulatedGrid::Private::Indices() const {
  return indices_;
}

const geom::Point3& TriangulatedGrid::Private::Origin() const {
  return origin_;
}

void TriangulatedGrid::Private::ApplyPose(const gz::math::Pose3d& pose) {
  // Origin - slide the patch in the xy - plane only
  geom::Point3 o = geom::Origin();
  origin_ = geom::Point3(o.x() + pose.Pos().X(), o.y() + pose.Pos().Y(), o.z());

  // Mesh points
  for (
    auto&& it = std::make_pair(std::begin(points0_), std::begin(points_));
    it.first != std::end(points0_) && it.second != std::end(points_);
    ++it.first, ++it.second)
  {
    const auto& p0 = *it.first;
    auto& p = *it.second;
    p = geom::Point3(p0.x() + pose.Pos().X(), p0.y() + pose.Pos().Y(), p0.z());
  }
}

bool TriangulatedGrid::Private::IsValid(bool verbose) const {
  bool isValid = true;
  const int64_t nPoints = static_cast<int64_t>(points_.size());
  if (nPoints != (nx_ + 1) * (ny_ + 1)) {
    isValid = false;
    if (verbose)
      std::cerr << "TriangulatedGrid: expected " << (nx_ + 1) * (ny_ + 1)
          << " points, have " << nPoints << "\n";
  }
  if (static_cast<int64_t>(indices_.size()) != 2 * nx_ * ny_) {
    isValid = false;
    if (verbose)
      std::cerr << "TriangulatedGrid: expected " << 2 * nx_ * ny_
          << " faces, have " << indices_.size() << "\n";
  }
  for (const auto& f : indices_) {
    for (int k = 0; k < 3; ++k) {
      if (f[k] < 0 || f[k] >= nPoints) {
        isValid = false;
        if (verbose)
          std::cerr << "TriangulatedGrid: face index out of range: "
              << f[k] << "\n";
      }
    }
  }
  return isValid;
}

void TriangulatedGrid::Private::DebugPrintMesh() const {
  std::cout << "nx: " << nx_ << "\n";
  std::cout << "ny: " << ny_ << "\n";
  std::cout << "lx: " << lx_ << "\n";
  std::cout << "ly: " << ly_ << "\n";

  std::cout << "points: [" << points_.size() << "]" << "\n";
  for (auto&& p : points_) {
    std::cout << p << "\n";
  }
  std::cout << "indices: [" << indices_.size() << "]" << "\n";
  for (auto&& i : indices_) {
    std::cout << i[0] << " " << i[1] << " " << i[2]  << "\n";
  }
  std::cout << "infinite indices: [" << infinite_indices_.size()
      << "]" << "\n";
  for (auto&& i : infinite_indices_) {
    std::cout << i[0] << " " << i[1] << "\n";
  }
}

void TriangulatedGrid::Private::DebugPrintTriangulation() const {
  std::cout << "regular lattice triangulation" << "\n";
  std::cout << "is valid: " << IsValid() << "\n";
  std::cout << "dimension: 2" << "\n";
  std::cout << "number of vertices : " << points_.size() << "\n";
  std::cout << "number of faces : " << indices_.size() << "\n";
  std::cout << "number of boundary edges : " << infinite_indices_.size()
      << "\n";
  std::cout << "origin : " << origin_ << "\n";
}

void TriangulatedGrid::Private::UpdatePoints(
      const std::vector<geom::Point3>& from) {
  points_ = from;
}

void TriangulatedGrid::Private::UpdatePoints(
    const std::vector<gz::math::Vector3d>& from) {
  auto it_to = points_.begin();
  auto it_from = from.begin();
  for ( ; it_to != points_.end() && it_from != from.end();
      ++it_to, ++it_from) {
    *it_to = geom::Point3(it_from->X(), it_from->Y(), it_from->Z());
  }
}

void TriangulatedGrid::Private::UpdatePoints(const geom::Mesh& from) {
  const Index n = std::min<Index>(points_.size(), geom::VertexCount(from));
  for (Index i = 0; i < n; ++i) {
    const geom::Point3& p = geom::VertexPoint(from, i);
    points_[i] = geom::Point3(p.x(), p.y(), p.z());
  }
}

TriangulatedGrid::~TriangulatedGrid()
{
}

//////////////////////////////////////////////////
TriangulatedGrid::TriangulatedGrid(Index nx, Index ny, double lx, double ly) :
  impl_(new TriangulatedGrid::Private(nx, ny, lx, ly))
{
}

//////////////////////////////////////////////////
void TriangulatedGrid::CreateMesh()
{
  impl_->CreateMesh();
}

//////////////////////////////////////////////////
void TriangulatedGrid::CreateTriangulation()
{
  impl_->CreateTriangulation();
}

//////////////////////////////////////////////////
std::unique_ptr<TriangulatedGrid> TriangulatedGrid::Create(
    Index nx, Index ny, double lx, double ly)
{
  std::unique_ptr<TriangulatedGrid> instance =
      std::make_unique<TriangulatedGrid>(nx, ny, lx, ly);
  instance->CreateMesh();
  instance->CreateTriangulation();
  return instance;
}

//////////////////////////////////////////////////
bool TriangulatedGrid::Locate(const geom::Point3& query,
    int64_t& faceIndex) const
{
  return impl_->Locate(query, faceIndex);
}

//////////////////////////////////////////////////
bool TriangulatedGrid::Height(const geom::Point3& query,
    double& height) const
{
  return impl_->Height(query, height);
}

//////////////////////////////////////////////////
bool TriangulatedGrid::Height(const std::vector<geom::Point3>& queries,
    std::vector<double>& heights) const
{
  return impl_->Height(queries, heights);
}

//////////////////////////////////////////////////
bool TriangulatedGrid::Interpolate(TriangulatedGrid& patch) const
{
  return impl_->Interpolate(patch);
}

//////////////////////////////////////////////////
const Point3Range& TriangulatedGrid::Points() const
{
  return impl_->Points();
}

//////////////////////////////////////////////////
const Index3Range& TriangulatedGrid::Indices() const
{
  return impl_->Indices();
}

//////////////////////////////////////////////////
const geom::Point3& TriangulatedGrid::Origin() const
{
  return impl_->Origin();
}

//////////////////////////////////////////////////
void TriangulatedGrid::ApplyPose(const gz::math::Pose3d& pose)
{
  impl_->ApplyPose(pose);
}

//////////////////////////////////////////////////
bool TriangulatedGrid::IsValid(bool verbose) const
{
  return impl_->IsValid(verbose);
}

//////////////////////////////////////////////////
void TriangulatedGrid::DebugPrintMesh() const
{
  impl_->DebugPrintMesh();
}

//////////////////////////////////////////////////
void TriangulatedGrid::DebugPrintTriangulation() const
{
  impl_->DebugPrintTriangulation();
}

//////////////////////////////////////////////////
void TriangulatedGrid::UpdatePoints(const std::vector<geom::Point3>& from)
{
  impl_->UpdatePoints(from);
}

//////////////////////////////////////////////////
void TriangulatedGrid::UpdatePoints(
    const std::vector<gz::math::Vector3d>& from)
{
  impl_->UpdatePoints(from);
}

//////////////////////////////////////////////////
void TriangulatedGrid::UpdatePoints(const geom::Mesh& from)
{
  impl_->UpdatePoints(from);
}

//////////////////////////////////////////////////
std::array<double, 2> TriangulatedGrid::TileSize() const
{
  return {impl_->lx_, impl_->ly_};
}

//////////////////////////////////////////////////
std::array<Index, 2> TriangulatedGrid::CellCount() const
{
  return {impl_->nx_, impl_->ny_};
}

}  // namespace waves
}  // namespace gz
