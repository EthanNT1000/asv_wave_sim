# CGAL usage audit

Scope: every `#include <CGAL/...>` and `CGAL::` symbol in `gz-waves/` (library,
Gazebo systems, unit tests). Read-only audit; no code was changed.

License facts were taken from the `SPDX-License-Identifier` line of each
header in the installed CGAL 6.2.1 (the same package split applies to CGAL 5.x,
which Ubuntu Jammy CI installs) and cross-checked against the package overview
at <https://doc.cgal.org/latest/Manual/packages.html>. Every CGAL package is
dual-licensed "GPL-3.0-or-later OR LicenseRef-Commercial" or
"LGPL-3.0-or-later OR LicenseRef-Commercial"; below "GPL" / "LGPL" means the
open-source arm of that pair.

## Verdict up front

**The production usage is NOT LGPL-only.** The library that ships in
`libgz-waves1.so` (and therefore every system plugin that links it) depends on
three GPL packages:

| GPL package | Where it is load-bearing | Can it be removed? |
|---|---|---|
| **Surface Mesh** (`CGAL::Surface_mesh`) | container for every hull collision mesh (`MeshTools::MakeSurfaceMesh`, `Physics.cc`, `Hydrodynamics.cc`) and for the water patch (`Grid`) | yes, hand-rollable: only "indexed triangle soup" features are used |
| **2D Triangulations** (`Constrained_Delaunay_triangulation_2`, `Triangulation_hierarchy_2`) | `TriangulatedGrid` point location behind `Wavefield::Height`, called every physics step from `WavefieldSampler::UpdatePatch` | yes, replaceable by index math on the regular ocean tile |
| **3D Fast Intersection and Distance Computation** (`CGAL::AABB_tree`) | `Geometry::MakeAABBTree` / `Geometry::SearchMesh`, compiled into the library but **no production caller** (tests only) | yes, dead code in production |

Everything else the code touches (kernel, number utilities, BGL helpers,
Timer, generators, STL extensions) is LGPL. So the CGAL footprint splits into
"GPL: three containers/structures that a regular-grid + triangle-soup
implementation does not need" and "LGPL: plain double-precision vector
arithmetic that a 30-line header could replace".

This matters only if the project's own GPL-3.0 licence is ever to change
(`LICENSE` is GPL-3.0; `README.md` already states CGAL and FFTW are GPL).
Note also that `LICENSE_THIRDPARTY` lists the Kernel (LGPL), Surface Mesh
(GPL) and AABB tree (GPL) but **omits the 2D Triangulations package (GPL)**
that `TriangulatedGrid.cc` uses; that attribution gap is a documentation fix
for a later change, not made here.

---

## 1. Inventory

Column "Use" abbreviations: **prod** = compiled into `libgz-waves1` and
reachable at runtime; **prod (dead)** = compiled in but no runtime caller;
**test** = `*_TEST.cc` only.

### 1a. GPL-dependent usage

| Include / symbol | CGAL package | License | Files / functions | What it does geometrically | Use |
|---|---|---|---|---|---|
| `<CGAL/Surface_mesh.h>`, `CGAL::Surface_mesh<Point3>` (`cgal::Mesh`), `Mesh::Vertex_index`, `Face_index`, `Halfedge_index`, `add_vertex`, `add_face`, `point()`, `vertices()`, `faces()`, `halfedge()/next()/target()`, `num_vertices/num_faces`, `add_property_map` | Surface Mesh | **GPL** | `CGALTypes.hh` typedefs; `MeshTools::MakeSurfaceMesh` (gz mesh → Surface_mesh); `Grid::Grid`, `Grid::GetPoint/SetPoint/GetTriangle/GetFace/RecalculateNormals`, `GridTools::FindIntersectionTriangle`; `Geometry::MakeTriangle`, `Geometry::Normal(mesh, face)`; `Physics.cc`: `Hydrodynamics::Hydrodynamics` (per-vertex depth property map `v:depth`), `UpdateSubmergedTriangles` (face loop, `vertices_around_face`); `Hydrodynamics.cc`: `ApplyPose` (rigid transform of all vertices), `CreateAxisAlignedBox`, mesh deep copy `std::make_shared<cgal::Mesh>(*initLinkMesh)`; `WavefieldSampler::ApplyPose/UpdatePatch` (patch vertex loop); `TriangulatedGrid::UpdatePoints(const Mesh&)` | Half-edge polygon mesh used purely as a vertex array + face→3 vertex lookup (triangle soup). No topological operation (no edge flips, no adjacency queries, no deletions) is performed anywhere. | prod |
| `<CGAL/Constrained_Delaunay_triangulation_2.h>`, `<CGAL/Constrained_triangulation_2.h>`, `<CGAL/Triangulation_2.h>`, `<CGAL/Triangulation_hierarchy_2.h>`, `<CGAL/Triangulation_face_base_with_info_2.h>`, `<CGAL/Triangulation_vertex_base_with_info_2.h>`; `CGAL::Constrained_Delaunay_triangulation_2`, `Triangulation_hierarchy_2`, `Triangulation_hierarchy_vertex_base_2`, `Triangulation_vertex_base_with_info_2`, `Constrained_triangulation_face_base_2`, `Triangulation_face_base_with_info_2`, `Triangulation_data_structure_2`, `No_constraint_intersection_tag`; methods `insert`, `insert_constraints`, `locate`, `set_point`, `finite_vertices_begin/end`, `all_faces_begin/end`, `infinite_vertex`, `is_valid`, `tds()` | 2D Triangulations (+ 2D Triangulation Data Structure) | **GPL** | `TriangulatedGrid::Private` (`CreateTriangulation`, `Locate`, `Height` ×2, `Interpolate`, `ApplyPose`, `UpdatePoints` ×3, `IsValid`, `DebugPrintTriangulation`); reached in production via `Wavefield::Height` ← `WavefieldSampler::UpdatePatch` (every step, 25 patch vertices per link) | Builds a constrained Delaunay triangulation of the `(nx+1)(ny+1)` ocean-tile lattice projected on the xy-plane, with the cell diagonals as constraints, plus a Kirkpatrick-style location hierarchy; `locate(p)` walks to the face containing the xy of the query, then a ray/triangle intersection gives the surface height. Vertex positions are later overwritten (`set_point`) with the displaced wave positions **without re-triangulating**. | prod |
| `<CGAL/Regular_triangulation_2.h>`, `CGAL::Regular_triangulation_2` | 2D Triangulations | **GPL** | `TriangulatedGrid.cc` (include only, unused), `CGAL_TEST.cc` | weighted Delaunay triangulation | test / dead include |
| `<CGAL/AABB_tree.h>`, `<CGAL/AABB_traits_3.h>` (`AABB_traits.h` for CGAL < 6), `<CGAL/AABB_face_graph_triangle_primitive.h>`; `CGAL::AABB_tree`, `AABB_traits_3`, `AABB_face_graph_triangle_primitive`, `first_intersection`, `Intersection_and_primitive_id` | 3D Fast Intersection and Distance Computation (AABB tree) | **GPL** | `CGALTypes.hh` (`cgal::AABBTree`); `Geometry::MakeAABBTree`, `Geometry::SearchMesh` ×2 (`Geometry.cc`); callers: `Geometry_TEST.cc`, `CGAL_TEST.cc` only | Bounding-volume hierarchy over mesh faces; first ray/mesh intersection in both ray directions. | prod (dead) / test |
| `<CGAL/Polyhedron_3.h>`, `CGAL::Polyhedron_3`, `CGAL::make_tetrahedron` (the latter from BGL generators, LGPL) | 3D Polyhedral Surface | **GPL** | `CGAL_TEST.cc` (`AABBPolyhedronFacetIntersection`) | half-edge polyhedron used as AABB primitive source | test |

### 1b. LGPL-only usage

| Include / symbol | CGAL package | License | Files / functions | What it does geometrically | Use |
|---|---|---|---|---|---|
| `<CGAL/Simple_cartesian.h>`, `CGAL::Simple_cartesian<double>` and its types `Point_3`, `Vector_3`, `Vector_2`, `Direction_2/3`, `Line_3`, `Ray_3`, `Triangle_3` (`cgal::Point3`, `Vector3`, `Vector2`, `Direction3`, `Line`, `Ray`, `Triangle`) | 2D and 3D Linear Geometry Kernel (Kernel_23) | LGPL | `CGALTypes.hh`; all of `Physics.hh/.cc`, `Geometry.cc`, `Grid.cc`, `WavefieldSampler.cc`, `WaterCurrentGrid.hh/.cc`, `Convert.cc`, `OceanTile.cc` (`OceanTilePrivate<cgal::Point3>`), `Wavefield.cc`, `Hydrodynamics.cc` | Plain double-precision 2D/3D points, vectors, directions, lines, rays and triangles with operator arithmetic; `Triangle_3` is a 3-point container plus `operator[]`. | prod |
| `CGAL::ORIGIN`, `CGAL::NULL_VECTOR` | Kernel_23 (`Origin.h`) | LGPL | 30 + 33 sites in `Physics.hh/.cc`, `Geometry.cc`, `Grid.cc`, `WavefieldSampler.cc`, `TriangulatedGrid.cc`, `Hydrodynamics.cc` | zero point / zero vector constants and point−origin conversion | prod |
| `CGAL::cross_product`, `CGAL::scalar_product` | Kernel_23 global functions | LGPL | `Geometry::TriangleArea`, `Ray/LineIntersectsTriangle`; `Physics.cc` (`ComputePointVelocities`, buoyancy/drag/lift torques `xr × F`, `cosTheta`, orientation checks in `Split*`/`AddFullySubmergedTriangle`); `Hydrodynamics.cc` aero pass | 3D cross and dot products | prod |
| `CGAL::normal(p,q,r)` | Kernel_23 (`global_functions_3`) | LGPL | `Geometry::Normal` ×3 | unnormalised triangle normal `(q−p)×(r−p)` | prod |
| `CGAL::collinear(p,q,r)` | Kernel_23 predicate | LGPL | `Geometry::Normal(const Triangle&)` | exact-zero test of the 3D orientation of three points (inexact with this kernel) | prod |
| `CGAL::centroid(p,q,r)`, `CGAL::midpoint(p,q)` | Kernel_23 (`global_functions`) | LGPL | `Geometry::TriangleCentroid` ×2, `Geometry::MidPoint` | `(p+q+r)/3`, `(p+q)/2` | prod |
| `CGAL::to_double` (`<CGAL/number_utils.h>`) | Number Types | LGPL | 26 sites in `Physics.cc`, `Grid.cc`, `Hydrodynamics.cc` | identity conversion for `double` field type | prod |
| `CGAL::Projection_traits_xy_3<Kernel>` (`<CGAL/Projection_traits_xy_3.h>`) | Kernel_23 (projection traits) | LGPL | `TriangulatedGrid::Private::Gt` (traits of the GPL triangulation above) | makes 3D points behave as their xy projection for 2D predicates | prod (only meaningful with the GPL triangulation) |
| `<CGAL/Exact_predicates_inexact_constructions_kernel.h>` | Kernel_23 / Filtered kernel | LGPL | `TriangulatedGrid.cc` include; the typedef using it is commented out; `CGAL_TEST.cc` | filtered exact predicates over doubles | dead include / test |
| `CGAL::vertices_around_face` (`Surface_mesh.h` pulls in `<CGAL/boost/graph/iterator.h>`) | CGAL and the Boost Graph Library (BGL) | LGPL | `Hydrodynamics::UpdateSubmergedTriangles` (`Physics.cc`) | iterate the three vertices of a face via half-edges (on the GPL container) | prod |
| `<CGAL/boost/graph/Euler_operations.h>`, `<CGAL/boost/graph/generators.h>`, `CGAL::make_tetrahedron` | BGL | LGPL | `CGAL_TEST.cc` | mesh construction helpers | test |
| `<CGAL/Timer.h>`, `CGAL::Timer` | Profiling tools, Hash Map, Union-find, Modifiers | LGPL | included in `Physics.cc`, `Geometry.cc`, `TriangulatedGrid.cc` (all uses commented out); live in `MeshTools_TEST.cc`, `CGAL_TEST.cc` | wall-clock timer | dead include / test |
| `<CGAL/algorithm.h>`, `CGAL::cpp11::copy_n` | STL Extensions | LGPL | `TriangulatedGrid.cc` (include only), `CGAL_TEST.cc` | `std::copy_n` shim | dead include / test |
| `<CGAL/point_generators_2.h>`, `CGAL::Random_points_in_square_2`, `CGAL::Creator_uniform_2` | Geometric Object Generators | LGPL | `TriangulatedGrid.cc` (include only), `CGAL_TEST.cc` | random 2D point generation for triangulation tests | dead include / test |
| `find_package(CGAL COMPONENTS Core)` | CGAL_Core (CORE number types) | LGPL | `gz-waves/CMakeLists.txt` | arbitrary-precision expression numbers; nothing in the code uses `CORE::`, the component is requested but unused | build only |

Not CGAL but pulled in by it: GMP/MPFR (LGPL-3.0 / GPL-2.0 dual) are linked
through `CGAL::CGAL` because the CGAL CMake target requires them; with
`Simple_cartesian<double>` no GMP arithmetic is ever executed.

---

## 2. Kernel

All production geometry uses one kernel, declared once:

```cpp
typedef CGAL::Simple_cartesian<double> Kernel;   // CGALTypes.hh
```

* **Inexact predicates and inexact constructions**, plain `double`, no filtering.
* `Exact_predicates_exact_constructions_kernel` is **not** used anywhere, so the
  question "does the code need exact constructions" is answered by the code
  itself: it never had them. All constructions (cross products, centroids,
  waterline intercepts, centres of pressure, ray/triangle parameters) are
  evaluated in floating point and are numerically insensitive (they are
  linear or bilinear in the inputs, no nested constructions are fed back into
  predicates except in the two cases listed in Sec. 5).
* `TriangulatedGrid.cc` includes `Exact_predicates_inexact_constructions_kernel.h`
  but the corresponding typedef is commented out; the triangulation traits are
  `Projection_traits_xy_3<Simple_cartesian<double>>`, i.e. **the CDT runs on
  inexact predicates**.

What switching the triangulation (or the whole code base) to
`Exact_predicates_inexact_constructions_kernel` (EPICK) would cost:

| Aspect | Cost |
|---|---|
| API | none: EPICK is `Filtered_kernel<Simple_cartesian<double>>`; `Point_3`, `Vector_3`, operators and all global functions have the same signatures, and `to_double` is already used at every extraction site. |
| Constructions | unchanged (still `double`). |
| Predicates (`collinear`, orientation/side-of-oriented-circle inside `locate`, `insert`, `insert_constraints`) | interval-arithmetic filter first, exact `Gmpq`/`MP_Float` fallback only when the filter cannot decide. Typically 1.5–3× the cost of the raw double predicate; in `locate` the predicate is the whole cost, so expect roughly 2× on `Wavefield::Height`. That path is called for 25 patch vertices per link per step, so the absolute cost is microseconds. |
| Force loop (`Physics.cc`) | no predicates are evaluated there except `CGAL::collinear` in `Geometry::Normal(Triangle)`, which is called once per face and sub-face; the rest is constructions, so no measurable change. |
| Compile time / binary | more template instantiation (filters, `Interval_nt`, `Lazy`); tens of seconds of compile time, no new link dependency beyond GMP/MPFR which are already linked. |
| Benefit | guaranteed-consistent orientation answers in the CDT walk (Sec. 5), correct detection of exactly degenerate triangles. |

Recommendation implied by Sec. 4: the only place exactness would buy anything is
the triangulation, and the triangulation is itself redundant.

---

## 3. Redundant usage given the regular grid

The ocean tile is an `(nx+1)×(ny+1)` lattice with spacing `lx/nx`, `ly/ny`
(`TriangulatedGrid::CreateMesh`, `Grid::Grid`, `OceanTile`). Each cell is split
by the same diagonal (`idx0-idx2`). For such a grid, point location is
`ix = floor((x - x0)/dx)`, `iy = floor((y - y0)/dy)`, and the sub-triangle is
decided by one comparison against the diagonal; interpolation is bilinear (or
barycentric on the chosen triangle, which is what the current ray/triangle
code computes anyway).

| Current CGAL-based path | What it computes | Replacement class | Notes |
|---|---|---|---|
| `TriangulatedGrid` CDT + hierarchy + `locate` → `Geometry::LineIntersectsTriangle` (`Wavefield::Height`) | wave-surface height at world `(x, y)` on the ocean tile | **replaceable by index math** | `GridTools::FindIntersectionIndex` already *is* this index math for the water patch; the tile path duplicates the lookup with a GPL structure. One caveat: for the `trochoid`/Gerstner and FFT algorithms the tile vertices are **horizontally displaced** (`OceanTile::UpdateVertex` adds `sx, sy`), and `TriangulatedGrid::UpdatePoints` moves the triangulation vertices to the displaced positions without re-triangulating, so the CDT is effectively used as an index into a mildly distorted lattice. Index math on the undisplaced lattice ignores that shift (error ≤ horizontal displacement, bounded by the steepness parameter); if the shift must be honoured, `WavefieldSampler::ComputeDepthDirectly` already contains the 2×2 Newton inversion of the Gerstner map, which is still index math plus a few iterations. Either way no triangulation is needed. |
| `Wavefield::Height` tile wrap-around (`fmod`) | periodic tiling | already index math | unchanged |
| `WavefieldSampler::ComputeDepth` → `GridTools::FindIntersectionIndex` + `FindIntersectionGrid` (shell search) + `LineIntersectsTriangle` (Möller–Trumbore on `cgal::Vector3`) | depth of a hull vertex below the 4×4 water patch | **already index math; CGAL is only the vector type → hand-rollable** | the expanding-shell fallback exists only because a ray/triangle test with an epsilon can miss a query exactly on a shared edge; bilinear/barycentric interpolation of the located cell is total and needs no fallback. |
| `Grid` (`Surface_mesh` as a regular-grid container; `std::advance` on vertex iterators to reach vertex `i`) | store patch vertices and face normals | **replaceable by index math** (`std::vector<Point>` of size `(nx+1)(ny+1)`, vertex `i` is `O(1)`) | `Grid::GetPoint/SetPoint/GetTriangle/GetFace` advance a bidirectional iterator `i` steps: `O(i)` per access. |
| `Hydrodynamics::PopulateSubmergedTriangle` / `SplitPartiallySubmergedTriangle1/2` / `AddFullySubmergedTriangle` (`Physics.cc`) | clip one hull triangle against the (locally flat) water surface using per-vertex heights | **already hand-rolled** (this *is* single-plane Sutherland–Hodgman on one triangle: 0, 1 or 2 output triangles); CGAL provides only `Point3`/`Vector3`/`Triangle` storage and `cross_product` for the orientation check → hand-rollable | nothing to remove except the types |
| `Physics::BuoyancyForceAtCenterOfPressure`, `CenterOfPressureApexUp/Dn`, `CenterOfForce`, `Geometry::HorizontalIntercept` | hydrostatic force and centre of pressure per clipped triangle | pure arithmetic → **hand-rollable** | uses only `+`, `-`, `*` on points/vectors |
| `Geometry::TriangleArea/Centroid/MidPoint/Normal/Normalize` | elementary triangle quantities | **hand-rollable** (`Kernel_23` global functions are one-liners) | |
| `Surface_mesh` as hull container (`MeshTools::MakeSurfaceMesh`, `ApplyPose`, `CreateAxisAlignedBox`, depth property map, `vertices_around_face`) | indexed triangle soup with a per-vertex scalar | **hand-rollable** (`struct { std::vector<Vec3> v; std::vector<std::array<int,3>> f; }`; `gz::common::Mesh` already holds exactly this and is the source of every `Surface_mesh` built) | the code relies on "vertex/face indices are 0..N−1 and stable", i.e. it already treats the mesh as arrays |
| `Geometry::MakeAABBTree` / `SearchMesh` (AABB tree) | first ray/mesh hit | **dead in production**; if revived: **needs a library** or a small BVH (a flat BVH over ~400–800 hull faces is ~150 lines; brute force over 800 faces is also acceptable at a few queries per step) | only `Geometry_TEST`/`CGAL_TEST` call it |
| `CGAL::Timer` | profiling | **hand-rollable** (`std::chrono::steady_clock`) | commented out in production |
| `CGAL::to_double` | no-op on `double` | delete | |

Summary: with a regular grid and triangle-soup hulls, **no CGAL data structure
in this repo is doing work that index arithmetic or a few lines of vector
algebra cannot do**. The single genuinely non-trivial algorithm (ray/mesh BVH)
has no production caller.

---

## 4. Where exact predicates are genuinely required

Ordered by how robustness-critical the path is.

1. **`Triangulation_2::locate` and CDT construction in `TriangulatedGrid`
   (production, every step).** The lattice is maximally degenerate for a
   Delaunay triangulation: the four corners of every cell are cocircular, and
   grid points are collinear along rows, columns and diagonals. `insert` on
   cocircular points and the `locate` walk both rely on `orientation` and
   `side_of_oriented_circle` predicates; with `Simple_cartesian<double>` these
   are evaluated in floating point. Two things keep it working today: lattice
   coordinates are exactly representable (`ix*dlx + lxm` with small integers),
   and the diagonal constraints pin the triangulation so no Delaunay flip is
   free to be decided by a rounding error. Both guarantees are lost once
   `UpdatePoints` overwrites the vertices with wave-displaced positions
   (Gerstner/FFT `sx, sy`): the structure is no longer Delaunay, may contain
   inverted faces for steep waves, and `locate` walks over it with inexact
   orientation tests. CGAL documents that an inexact kernel can make `locate`
   return a wrong face or fail to terminate on degenerate input. **If the
   triangulation stays, this is the one place that genuinely needs EPICK (and
   a re-triangulation after displacement); the cheaper fix is Sec. 3, drop the
   triangulation.**
2. **`CGAL::collinear` in `Geometry::Normal(const Triangle&)`.** Used to return a
   `NULL_VECTOR` for degenerate triangles. Degenerate sub-triangles are
   *common* here: `SplitPartiallySubmergedTriangle1/2` produce a zero-area
   triangle whenever a vertex lies exactly at the surface (`hl == 0`), and
   `HorizontalIntercept` produces `D == L` for horizontal faces. With an inexact
   kernel `collinear` is an exact-zero test on a double determinant, so it
   catches exactly-degenerate inputs but not near-degenerate slivers; the
   callers guard with `n == NULL_VECTOR` and `squared_length` checks, so the
   consequence is a slightly noisy normal on a face of ~zero area, i.e. ~zero
   force. Exact predicates would make the classification consistent; a
   tolerance on the area is the practical fix. Robustness-relevant, not
   critical.
3. **Ray/triangle tests on shared edges (`LineIntersectsTriangle` in
   `GridTools::FindIntersectionTriangle` and `TriangulatedGrid::Height`).**
   A query whose xy lies exactly on a cell diagonal or edge can fail the
   `u, v ∈ [0,1]` test in both adjacent triangles after rounding. The code
   mitigates with `FindIntersectionCell` (try both triangles) and the
   expanding-shell search; an exact predicate would remove the miss, an
   inclusive barycentric interpolation removes the problem entirely.
4. **Waterline classification `hh > 0`, `hm > 0`, `hl > 0`** in
   `PopulateSubmergedTriangle`. Not a geometric predicate: each vertex height is
   a single double compared with zero, and the same value is used for every
   face sharing the vertex, so adjacent faces classify consistently. No exact
   arithmetic needed.
5. **`AABB_tree::first_intersection`** (dead in production). Rays hitting a
   shared edge between two hull faces can be reported twice or missed with an
   inexact kernel; irrelevant while the function has no runtime caller.

Nothing in the force computation (`Physics.cc` buoyancy, drag, lift, damping)
requires exact predicates: it consists of constructions whose errors are
bounded by machine epsilon times the magnitudes involved, and the results are
summed forces, not combinatorial decisions.

---

## 5. Test-only CGAL usage (for completeness)

`CGAL_TEST.cc` (15 tests) exercises `Polyhedron_3` + AABB tree, `Surface_mesh`
grids, plain/constrained/hierarchy triangulations, `Regular_triangulation_2`,
random point generators and `Timer`; `TriangulatedGrid_TEST.cc` (3),
`Grid_TEST.cc` (5), `Geometry_TEST.cc` (8) test the wrappers above. These are
compiled only with `BUILD_TESTING=ON` and never ship in the library. They are
also the only exercisers of `Polyhedron_3` and `Regular_triangulation_2`.

---

## Appendix: file → CGAL dependency map

| File | GPL packages | LGPL packages |
|---|---|---|
| `include/gz/waves/CGALTypes.hh` | Surface Mesh, AABB tree | Kernel_23 |
| `src/Physics.cc`, `include/.../Physics.hh` | Surface Mesh (container, property map) | Kernel_23, Number Types, BGL (`vertices_around_face`), Timer (unused include) |
| `src/Geometry.cc` | Surface Mesh (face → points), AABB tree (dead) | Kernel_23, Timer (unused include) |
| `src/Grid.cc` | Surface Mesh | Kernel_23, Number Types |
| `src/TriangulatedGrid.cc` | 2D Triangulations, (Regular triangulation include unused) | Kernel_23 (projection traits, EPICK include unused), STL Extensions, Generators, Timer (unused includes) |
| `src/WavefieldSampler.cc` | Surface Mesh (via `Grid`) | Kernel_23 |
| `src/Wavefield.cc` | 2D Triangulations (via `TriangulatedGrid`) | Kernel_23 |
| `src/MeshTools.cc` | Surface Mesh | Kernel_23 |
| `src/OceanTile.cc`, `src/Convert.cc`, `src/WaterCurrentGrid.cc` | – | Kernel_23 (point/vector types only) |
| `src/systems/hydrodynamics/Hydrodynamics.cc` | Surface Mesh (copy, transform, bbox) | Kernel_23, Number Types |
| `src/CGAL_TEST.cc` | Surface Mesh, AABB tree, Polyhedron, 2D Triangulations | Kernel_23, BGL, Generators, STL Extensions, Timer |
