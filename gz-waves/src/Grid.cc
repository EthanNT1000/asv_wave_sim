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

#include "gz/waves/Grid.hh"

#include <algorithm>
#include <array>
#include <memory>
#include <iostream>
#include <iterator>
#include <string>
#include <vector>

#include <gz/common/Console.hh>

#include "gz/waves/geom/Geom.hh"
#include "gz/waves/Geometry.hh"
#include "gz/waves/Types.hh"

namespace gz
{
namespace waves
{

//////////////////////////////////////////////////
class GridPrivate
{
 public:
  /// \brief The size of the grid in each direction
  std::array<double, 2> size;

  /// \brief The number of cells in each direction
  std::array<Index, 2> cellCount;

  /// \brief The position of the grid center
  geom::Point3 center;

  /// \brief The grid mesh
  std::shared_ptr<geom::Mesh> mesh;

  /// \brief The grid normals (for each face)
  std::vector<geom::Vector3> normals;
};

//////////////////////////////////////////////////
Grid::Grid(
  const std::array<double, 2>& _size,
  const std::array<Index, 2>& _cellCount) :
    data(new GridPrivate())
{
  this->data->size = _size;
  this->data->cellCount = _cellCount;
  this->data->center = geom::Origin();
  this->data->mesh = std::make_shared<geom::Mesh>();

  // Grid dimensions
  const Index nx = this->data->cellCount[0];
  const Index ny = this->data->cellCount[1];
  const double Lx = this->data->size[0];
  const double Ly = this->data->size[1];
  const double lx = Lx/nx;
  const double ly = Ly/ny;

  auto& mesh = *this->data->mesh;
  auto& normals = this->data->normals;

  // Add vertices
  for (Index iy=0; iy <= ny; ++iy)
  {
    double py = iy * ly - Ly/2.0;
    for (Index ix=0; ix <= nx; ++ix)
    {
      double px = ix * lx - Lx/2.0;
      geom::AddVertex(mesh, geom::Point3(px, py, 0));
    }
  }

  // Add faces
  for (Index iy=0; iy < ny; ++iy)
  {
    for (Index ix=0; ix < nx; ++ix)
    {
      // Get the vertices in the cell coordinates
      const Index idx0 = iy * (nx+1) + ix;
      const Index idx1 = iy * (nx+1) + ix + 1;
      const Index idx2 = (iy+1) * (nx+1) + ix + 1;
      const Index idx3 = (iy+1) * (nx+1) + ix;

      // Faces
      geom::AddFace(mesh, idx0, idx1, idx2);
      geom::AddFace(mesh, idx0, idx2, idx3);

      // Face Normals
      geom::Point3 p0(geom::VertexPoint(mesh, idx0));
      geom::Point3 p1(geom::VertexPoint(mesh, idx1));
      geom::Point3 p2(geom::VertexPoint(mesh, idx2));
      geom::Point3 p3(geom::VertexPoint(mesh, idx3));
      normals.push_back(Geometry::Normal(p0, p1, p2));
      normals.push_back(Geometry::Normal(p0, p2, p3));
    }
  }
}

//////////////////////////////////////////////////
Grid::Grid(const Grid& _other) :
  data(new GridPrivate())
{
  this->data->size = _other.data->size;
  this->data->cellCount = _other.data->cellCount;
  this->data->mesh.reset(new geom::Mesh(*_other.data->mesh));
}

//////////////////////////////////////////////////
Grid& Grid::operator=(const Grid& _other)
{
  // Check for self assignment
  if (&_other == this)
    return *this;

  // Copy
  this->data->size = _other.data->size;
  this->data->cellCount = _other.data->cellCount;
  this->data->mesh.reset(new geom::Mesh(*_other.data->mesh));
  return *this;
}

//////////////////////////////////////////////////
std::shared_ptr<const geom::Mesh> Grid::GetMesh() const
{
  return this->data->mesh;
}

//////////////////////////////////////////////////
std::shared_ptr<geom::Mesh> Grid::GetMesh()
{
  return this->data->mesh;
}

//////////////////////////////////////////////////
const geom::Mesh& Grid::GetMeshByRef() const
{
  return *this->data->mesh;
}

//////////////////////////////////////////////////
const std::array<double, 2>& Grid::GetSize() const
{
  return this->data->size;
}

//////////////////////////////////////////////////
const std::array<Index, 2>& Grid::GetCellCount() const
{
  return this->data->cellCount;
}

//////////////////////////////////////////////////
Index Grid::GetVertexCount() const
{
  return geom::VertexCount(*this->data->mesh);
}

//////////////////////////////////////////////////
Index Grid::GetFaceCount() const
{
  return geom::FaceCount(*this->data->mesh);
}

//////////////////////////////////////////////////
const geom::Point3& Grid::GetPoint(Index _i) const
{
  return geom::VertexPoint(*this->data->mesh, _i);
}

//////////////////////////////////////////////////
void Grid::SetPoint(Index _i, const geom::Point3& _v)
{
  geom::SetVertexPoint(*this->data->mesh, _i, _v);
}

//////////////////////////////////////////////////
geom::Triangle Grid::GetTriangle(Index _ix, Index _iy, Index _k) const
{
  // Original lookup using cell indexing - keep for index arithmetic
  // // Grid dimensions
  // const Index nx = this->data->cellCount[0];
  // const Index ny = this->data->cellCount[1];

  // // Get the vertices in the cell coordinates
  // const Index idx0 = _iy * (nx+1) + _ix;
  // const Index idx1 = _iy * (nx+1) + _ix + 1;
  // const Index idx2 = (_iy+1) * (nx+1) + _ix + 1;
  // const Index idx3 = (_iy+1) * (nx+1) + _ix;

  // switch (_k)
  // {
  // case 0:
  //   return geom::Triangle(
  //     this->GetPoint(idx0),
  //     this->GetPoint(idx1),
  //     this->GetPoint(idx2)
  //   );
  // case 1:
  //   return geom::Triangle(
  //     this->GetPoint(idx0),
  //     this->GetPoint(idx2),
  //     this->GetPoint(idx3)
  //   );
  // default:
  //   gzthrow("Index too large: " << _k << " > 1");
  // }

  // Optimised lookup using face indexing
  // Face index
  const Index nx = this->data->cellCount[0];
  const Index idx = 2 * (nx * _iy + _ix) + _k;

  // Make triangle from face index
  return geom::FaceTriangle(*this->data->mesh, idx);
}

//////////////////////////////////////////////////
geom::FaceIndex Grid::GetFace(Index _ix, Index _iy, Index _k) const
{
  // Face index
  const Index nx = this->data->cellCount[0];
  const Index idx = 2 * (nx * _iy + _ix) + _k;

  return geom::ToFaceIndex(idx);
}

//////////////////////////////////////////////////
const geom::Vector3& Grid::GetNormal(Index _ix, Index _iy, Index _k) const
{
  // Face index
  const Index nx = this->data->cellCount[0];
  const Index idx = 2 * (nx * _iy + _ix) + _k;

  return this->data->normals[idx];
}

//////////////////////////////////////////////////
const geom::Vector3& Grid::GetNormal(Index _idx) const
{
  return this->data->normals[_idx];
}

//////////////////////////////////////////////////
void Grid::RecalculateNormals()
{
  auto& mesh = *this->data->mesh;
  const Index nFaces = geom::FaceCount(mesh);
  for (Index idx = 0; idx < nFaces; ++idx)
  {
    this->data->normals[idx] = Geometry::Normal(mesh, geom::ToFaceIndex(idx));
  }
}

//////////////////////////////////////////////////
const geom::Point3& Grid::GetCenter() const
{
  return this->data->center;
}

//////////////////////////////////////////////////
void Grid::SetCenter(const geom::Point3& _center)
{
  this->data->center = _center;
}

//////////////////////////////////////////////////
void Grid::DebugPrint() const
{
  auto&& mesh = *this->data->mesh;

  gzmsg << "Center " << std::endl;
  gzmsg << "c0:  " << this->data->center << std::endl;

  gzmsg << "Vertices " << std::endl;
  const Index nVerts = geom::VertexCount(mesh);
  for (Index v = 0; v < nVerts; ++v)
  {
    gzmsg << v << ": " << geom::VertexPoint(mesh, v) << std::endl;
  }
  gzmsg << "Faces " << std::endl;
  const Index nFaces = geom::FaceCount(mesh);
  for (Index f = 0; f < nFaces; ++f)
  {
    geom::Triangle tri = geom::FaceTriangle(mesh, f);
    gzmsg << f << ": " << tri << std::endl;
  }
}

//////////////////////////////////////////////////
//////////////////////////////////////////////////
bool GridTools::FindIntersectionIndex(
  const Grid& _grid,
  double _x, double _y,
  std::array<Index, 3>& _index
)
{
  const Index nx = _grid.GetCellCount()[0];
  const Index ny = _grid.GetCellCount()[1];
  const double Lx = _grid.GetSize()[0];
  const double Ly = _grid.GetSize()[1];
  const double lx = Lx/nx;
  const double ly = Ly/ny;
  const double lowerX = -Lx/2.0 + _grid.GetCenter().x();
  const double upperX =  Lx/2.0 + _grid.GetCenter().x();
  const double lowerY = -Ly/2.0 + _grid.GetCenter().y();
  const double upperY =  Ly/2.0 + _grid.GetCenter().y();

  // Check x bounds
  if (_x < lowerX || _x > upperX)
    return false;

  // Check y bounds
  if (_y < lowerY || _y > upperY)
    return false;

  // Cell
  Index ix = static_cast<Index>(std::floor((_x - lowerX)/lx));
  Index iy = static_cast<Index>(std::floor((_y - lowerY)/ly));

  // Face / triangle
  double x0 = ix * lx + lowerX;
  double y0 = iy * ly + lowerY;
  double m = ly/lx;
  double c = y0 - m * x0;
  Index k = _y > (m * _x + c) ? 1 : 0;

  // @DEBUG_INFO
  // gzmsg << "ix: " << ix << std::endl;
  // gzmsg << "iy: " << iy << std::endl;
  // gzmsg << "x0: " << x0 << std::endl;
  // gzmsg << "y0: " << y0 << std::endl;
  // gzmsg << "m:  " << m  << std::endl;
  // gzmsg << "c:  " << c  << std::endl;
  // gzmsg << "k:  " << k  << std::endl;

  _index[0] = ix;
  _index[1] = iy;
  _index[2] = k;

  return true;
}

//////////////////////////////////////////////////
bool GridTools::FindIntersectionTriangle(
  const Grid& _grid,
  const geom::Point3& _origin,
  const geom::Direction3& _direction,
  const std::array<Index, 3>& _index,
  geom::Point3& _intersection)
{
  // FaceIndex version: the ByRef vs shared_ptr access makes
  // a difference (50% of this functions execution time!)
  const auto& mesh = _grid.GetMeshByRef();
  auto face = _grid.GetFace(_index[0], _index[1], _index[2]);
  const auto v = geom::FaceVertices(mesh, geom::ToIndex(face));
  const geom::Point3& p0 = geom::VertexPoint(mesh, v[0]);
  const geom::Point3& p1 = geom::VertexPoint(mesh, v[1]);
  const geom::Point3& p2 = geom::VertexPoint(mesh, v[2]);

  return Geometry::LineIntersectsTriangle(
    _origin, _direction, p0, p1, p2, _intersection);
}

//////////////////////////////////////////////////
bool GridTools::FindIntersectionCell(
  const Grid& _grid,
  const geom::Point3& _origin,
  const geom::Direction3& _direction,
  std::array<Index, 3>& _index,
  geom::Point3& _intersection
)
{
  // Search each of the two triangles comprising each cell
  bool isFound = FindIntersectionTriangle(
    _grid, _origin, _direction, _index, _intersection);
  if (isFound)
  {
    return true;
  }
  _index[2] = (_index[2] == 0) ? 1 : 0;
  return FindIntersectionTriangle(
    _grid, _origin, _direction, _index, _intersection);
}

//////////////////////////////////////////////////
// NOTE: This routine as written must use 'Index' rather than 'unsigned int'
// otherwise the index arithmetic will be incorrect.
bool GridTools::FindIntersectionGrid(
  const Grid& _grid,
  const geom::Point3& _origin,
  const geom::Direction3& _direction,
  std::array<Index, 3>& _index,
  geom::Point3& _point
)
{
  // The flag 'isDone' is true if there are no remaining cells to search.
  bool isDone = false;

  // index of the first cell searched.
  isDone = FindIntersectionCell(_grid, _origin, _direction, _index, _point);
  if (isDone)
  {
    return isDone;
  }

  // Grid dimensions
  Index nx = _grid.GetCellCount()[0];
  Index ny = _grid.GetCellCount()[1];

  // indexes for the grid boundaries
  Index kxmin = 0;
  Index kxmax = nx - 1;
  Index kymin = 0;
  Index kymax = ny - 1;

  // grid boundary indexes for the first cell.
  Index kx   = _index[0];
  Index ky   = _index[1];
  Index kxm0 = kx;
  Index kxp0 = kxm0;
  Index kym0 = ky;
  Index kyp0 = kym0;

  // loop searches the cells in an expanding shell about
  // the area already searched.
  while (!isDone)
  {
    isDone = true;

    // boundary of the next shell
    Index kxm = std::max(kxm0 - 1, kxmin);
    Index kxp = std::min(kxp0 + 1, kxmax);
    Index kym = std::max(kym0 - 1, kymin);
    Index kyp = std::min(kyp0 + 1, kymax);

    // row0 : lower limit offset for the column loop is either 0 or 1
    Index row0 = 0;
    if (kxm != kxm0)
    {
      isDone = false;
      row0 = 1;
      for (Index j=kym; j <= kyp; ++j)
      {
        _index[0] = kxm;
        _index[1] = j;
        isDone = FindIntersectionCell(_grid, _origin,
            _direction, _index, _point);
        if (isDone)
          return true;
      }
    }

    // row1 : upper limit offset for the column loop is either 0 or 1.
    Index row1 = 0;
    if (kxp != kxp0)
    {
      isDone = false;
      row1 = 1;
      for (Index j=kym; j <= kyp; ++j)
      {
        _index[0] = kxp;
        _index[1] = j;
        isDone = FindIntersectionCell(_grid, _origin,
            _direction, _index, _point);
        if (isDone)
          return true;
      }
    }

    // col0
    if (kym != kym0)
    {
      isDone = false;
      for (Index i=kxm+row0; i <= kxp-row1; ++i)
      {
        _index[0] = i;
        _index[1] = kym;
        isDone = FindIntersectionCell(_grid, _origin,
            _direction, _index, _point);
        if (isDone)
          return true;
      }
    }

    // col1
    if (kyp != kyp0)
    {
      isDone = false;
      for (Index i=kxm+row0; i <= kxp-row1; ++i)
      {
        _index[0] = i;
        _index[1] = kyp;
        isDone = FindIntersectionCell(_grid, _origin,
            _direction, _index, _point);
        if (isDone)
          return true;
      }
    }

    // update
    kxm0 = kxm;
    kxp0 = kxp;
    kym0 = kym;
    kyp0 = kyp;
  }

  // Fail to find an intersection after searching all cells
  return false;
}

}  // namespace waves
}  // namespace gz
