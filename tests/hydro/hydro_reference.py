"""Pure-numpy reference port of the gz-waves hydrodynamics force model.

This module mirrors, term by term, the algorithms in
``gz-waves/src/Physics.cc`` (Kerner-style per-triangle buoyancy, ITTC-57
viscous drag, pressure drag, foil lift, per-DOF Fossen damping) together with
the Fossen matrix constructions the audit measures them against
(``docs/hydrodynamics_audit.md``).

It is used by ``tests/hydro/test_properties.py`` and ``scripts/decay_test.py``.
It has no dependency other than numpy so it runs without a Gazebo install.

Conventions
-----------
* World frame is z-up, gravity is ``-9.81`` along z, as in
  ``PhysicalConstants.cc``.
* Mesh triangle normals are *outward*, which is what the buoyancy term in
  ``Physics.cc`` requires (F = rho*g*A*h*n with g negative).
* Body-fixed velocity vector nu = [u v w p q r] follows Fossen (SNAME).
"""
from __future__ import annotations

import math
from dataclasses import dataclass, field

import numpy as np

# PhysicalConstants.cc
GRAVITY = -9.81                # m/s^2 (z-up)
RHO_WATER = 1025.0             # kg/m^3
NU_WATER = 1.0533e-6           # m^2/s
RHO_AIR = 1.225                # kg/m^3 (Hydrodynamics.cc)

# HydrodynamicsParametersPrivate defaults (Physics.cc)
DEFAULT_PARAMS = dict(
    cDampU1=1.0e-6, cDampU2=1.0e-6,
    cDampV1=1.0e-3, cDampV2=1.0e-3,
    cDampW1=1.0e-3, cDampW2=1.0e-3,
    cDampP1=5.0e-3, cDampP2=5.0e-3,
    cDampQ1=5.0e-3, cDampQ2=5.0e-3,
    cDampN1=5.0e-4, cDampN2=5.0e-4,
    cPDrag1=1.0e2, cPDrag2=1.0e2, fPDrag=0.4,
    cSDrag1=1.0e2, cSDrag2=1.0e2, fSDrag=0.4,
    vRDrag=1.0,
    cLift1=1.0, alphaStall=0.26,
)
DEFAULT_PARAMS["cLMax"] = (DEFAULT_PARAMS["cLift1"] * 2.0 * math.pi
                           * math.sin(DEFAULT_PARAMS["alphaStall"]))
BOTTOM_THRESHOLD = 0.3


# ---------------------------------------------------------------------------
# Fossen matrix constructions (Ch. 3 rigid-body kinetics, Ch. 6 added mass)
# ---------------------------------------------------------------------------
def skew(v):
    """S(v) such that S(v) a = v x a (Fossen Def. 2.2)."""
    x, y, z = v
    return np.array([[0.0, -z, y],
                     [z, 0.0, -x],
                     [-y, x, 0.0]])


def rigid_body_mass_matrix(m, inertia_cg, r_g=(0.0, 0.0, 0.0)):
    """M_RB about the body origin (Fossen Sec. 3.3, Eq. 3.44).

    ``inertia_cg`` is the 3x3 inertia tensor about the CG expressed in body
    axes, ``r_g`` the CG position relative to the body origin.
    """
    r_g = np.asarray(r_g, dtype=float)
    S = skew(r_g)
    I_o = np.asarray(inertia_cg, dtype=float) - m * S @ S   # parallel axis
    M = np.zeros((6, 6))
    M[:3, :3] = m * np.eye(3)
    M[:3, 3:] = -m * S
    M[3:, :3] = m * S
    M[3:, 3:] = I_o
    return M


def coriolis_from_mass_matrix(M, nu):
    """C(nu) derived from a (symmetric) mass matrix by the velocity-independent
    parametrization of Fossen Theorem 3.2 (Eq. 3.60 for M_RB, applied to M_A
    in Sec. 6.3 to give C_A(nu_r)).  The result is skew-symmetric for any nu.
    """
    M = np.asarray(M, dtype=float)
    nu = np.asarray(nu, dtype=float)
    M11, M12 = M[:3, :3], M[:3, 3:]
    M21, M22 = M[3:, :3], M[3:, 3:]
    nu1, nu2 = nu[:3], nu[3:]
    C = np.zeros((6, 6))
    C[:3, 3:] = -skew(M11 @ nu1 + M12 @ nu2)
    C[3:, :3] = -skew(M11 @ nu1 + M12 @ nu2)
    C[3:, 3:] = -skew(M21 @ nu1 + M22 @ nu2)
    return C


def munk_moment(M_A, nu_r):
    """Pitch/yaw moment of C_A(nu_r) nu_r (the Munk moment, Fossen Sec. 6.3)."""
    return (coriolis_from_mass_matrix(M_A, nu_r) @ nu_r)[3:]


# ---------------------------------------------------------------------------
# Per-DOF damping exactly as ComputeDampingForce (Physics.cc)
# ---------------------------------------------------------------------------
_DOF_KEYS = ["U", "V", "W", "P", "Q", "N"]


def damping_matrix(nu_r, params=DEFAULT_PARAMS):
    """D(nu_r) = D_L + D_n(nu_r) as a 6x6 diagonal matrix (Fossen Eq. 6.53).

    Mirrors ``Hydrodynamics::ComputeDampingForce``: per DOF
    F_i = -(c1_i * v_i + c2_i * |v_i| * v_i), i.e. D_ii = c1_i + c2_i |v_i|.
    """
    d = []
    for i, key in enumerate(_DOF_KEYS):
        c1 = params[f"cDamp{key}1"]
        c2 = params[f"cDamp{key}2"]
        d.append(c1 + c2 * abs(nu_r[i]))
    return np.diag(d)


def damping_force_world(R_bw, v_world, omega_world, v_current_world,
                        params=DEFAULT_PARAMS):
    """World-frame (force, torque) from ComputeDampingForce.

    ``R_bw`` rotates body vectors into the world frame.  Only the linear
    velocity is taken relative to the (irrotational) current, matching
    Fossen's nu_c = [u_c v_c w_c 0 0 0] (Sec. 8.3 / Eq. 8.x irrotational
    current) and the implementation.
    """
    lin_b = R_bw.T @ (np.asarray(v_world) - np.asarray(v_current_world))
    ang_b = R_bw.T @ np.asarray(omega_world)
    nu_r = np.concatenate([lin_b, ang_b])
    tau = -(damping_matrix(nu_r, params) @ nu_r)
    return R_bw @ tau[:3], R_bw @ tau[3:]


# ---------------------------------------------------------------------------
# Per-triangle hydrostatics (Physics::BuoyancyForceAtCenterOfPressure)
# ---------------------------------------------------------------------------
def _center_of_force(fA, fB, A, B):
    div = fA + fB
    if abs(div) > 1e-16:
        t = fA / div
        return B + (A - B) * t
    return B


def _cop_apex_up(z0, H, M, B):
    alt = B - H
    h = H[2] - M[2]
    tc = 2.0 / 3.0
    div = 6.0 * z0 + 4.0 * h
    if abs(div) > 1e-16:
        tc = (4.0 * z0 + 3.0 * h) / div
    return H + alt * tc


def _cop_apex_dn(z0, L, M, B):
    alt = L - B
    h = M[2] - L[2]
    tc = 1.0 / 3.0
    div = 6.0 * z0 + 2.0 * h
    if abs(div) > 1e-16:
        tc = (2.0 * z0 + h) / div
    return B + alt * tc


def _horizontal_intercept(H, M, L):
    """Point D on edge L-H at the height of M (Geometry::HorizontalIntercept).

    For a horizontal triangle (H.z == L.z) the C++ code keeps t = 0 and
    returns L (at M's height), so the full triangle is treated as apex-up.
    """
    t = 0.0
    div = H[2] - L[2]
    if abs(div) > np.finfo(float).eps:
        t = (M[2] - L[2]) / div
    return np.array([L[0] + t * (H[0] - L[0]), L[1] + t * (H[1] - L[1]), M[2]])


def triangle_normal(tri):
    n = np.cross(tri[1] - tri[0], tri[2] - tri[0])
    ln = np.linalg.norm(n)
    return n / ln if ln > 0 else np.zeros(3)


def triangle_area(tri):
    return 0.5 * np.linalg.norm(np.cross(tri[1] - tri[0], tri[2] - tri[0]))


def buoyancy_force_at_cop(tri, depth_fn, rho=RHO_WATER, g=GRAVITY):
    """Port of Physics::BuoyancyForceAtCenterOfPressure(sampler, triangle).

    ``depth_fn(point)`` returns the depth h below the local free surface
    (positive when submerged), as WavefieldSampler::ComputeDepth does.
    Returns (center, force).
    """
    tri = np.asarray(tri, dtype=float)
    order = np.argsort(-tri[:, 2])          # sort by height, highest first
    H, M, L = tri[order[0]], tri[order[1]], tri[order[2]]
    C = tri.mean(axis=0)
    depthC = depth_fn(C)
    normal = triangle_normal(tri)

    D = _horizontal_intercept(H, M, L)
    B = 0.5 * (M + D)
    fU = fL = 0.0
    CpU = B.copy()
    CpL = B.copy()

    if H[2] >= M[2]:
        z0 = depthC - (H[2] - C[2])
        CpU = _cop_apex_up(z0, H, M, B)
        CU = (H + M + D) / 3.0
        hCU = depthC + (C[2] - CU[2])
        area = triangle_area(np.array([H, M, D]))
        fU = rho * g * area * hCU
    if M[2] > L[2]:
        z0 = depthC + (C[2] - M[2])
        CpL = _cop_apex_dn(z0, L, M, B)
        CL = (L + M + D) / 3.0
        hCL = depthC + (C[2] - CL[2])
        area = triangle_area(np.array([L, M, D]))
        fL = rho * g * area * hCL

    force = normal * (fU + fL)
    center = _center_of_force(fU, fL, CpU, CpL)
    return center, force


def split_submerged(tri, heights):
    """Port of PopulateSubmergedTriangle / SplitPartiallySubmergedTriangle*.

    ``heights`` are vertex heights above the water surface (negative depth).
    Returns a list of submerged sub-triangles with the original orientation.
    """
    tri = np.asarray(tri, dtype=float)
    heights = np.asarray(heights, dtype=float)
    n = triangle_normal(tri)
    idx = np.argsort(-heights)              # H, M, L by height
    hh, hm, hl = heights[idx]
    vh, vm, vl = tri[idx]

    def orient(t):
        return t if np.dot(n, triangle_normal(t)) >= 0.0 else t[[0, 2, 1]]

    if hh > 0:
        if hm > 0:
            if hl > 0:
                return []
            tm = -hl / (hm - hl)
            th = -hl / (hh - hl)
            vmi = vl + (vm - vl) * tm
            vhi = vl + (vh - vl) * th
            return [orient(np.array([vl, vmi, vhi]))]
        tm = -hm / (hh - hm)
        tl = -hl / (hh - hl)
        vmi = vm + (vh - vm) * tm
        vli = vl + (vh - vl) * tl
        return [orient(np.array([vm, vmi, vl])), orient(np.array([vmi, vli, vl]))]
    return [orient(np.array([vh, vm, vl]))]


def hydrostatic_wrench(tris, depth_fn, com, rho=RHO_WATER, g=GRAVITY):
    """Total buoyancy force and torque about ``com`` for a triangle soup.

    Also returns the waterline segments' bounding extents (for L_wl / beam)
    and the list of submerged sub-triangles, mirroring Hydrodynamics::Update.
    """
    force = np.zeros(3)
    torque = np.zeros(3)
    sub = []
    for tri in tris:
        tri = np.asarray(tri, dtype=float)
        heights = np.array([-depth_fn(v) for v in tri])
        for st in split_submerged(tri, heights):
            c, f = buoyancy_force_at_cop(st, depth_fn, rho, g)
            force += f
            torque += np.cross(c - com, f)
            sub.append(st)
    return force, torque, sub


# ---------------------------------------------------------------------------
# Per-triangle drag / lift terms (Physics.cc)
# ---------------------------------------------------------------------------
def ittc57_cf(Rn):
    r = max(1.0e3, Rn)
    d = math.log10(r) - 2.0
    return 0.075 / (d * d)


@dataclass
class TriKinematics:
    normal: np.ndarray
    centroid: np.ndarray
    area: float
    xr: np.ndarray
    v_rel: np.ndarray
    v_rel_mag: float
    up: np.ndarray
    cosTheta: float
    ut: np.ndarray
    alpha: float


def point_kinematics(tri, com, v_com, omega, v_fluid_fn):
    """Port of ComputePointVelocities using CoM velocity + omega x r."""
    tri = np.asarray(tri, dtype=float)
    n = triangle_normal(tri)
    c = tri.mean(axis=0)
    xr = c - com
    vp = np.asarray(v_com) + np.cross(omega, xr)
    v_rel = vp - v_fluid_fn(c)
    mag = np.linalg.norm(v_rel)
    vdn = float(np.dot(v_rel, n))
    v_rel_t = v_rel - n * vdn
    vt_mag = np.linalg.norm(v_rel_t)
    alpha = math.atan2(vdn, vt_mag + 1e-9)
    up = v_rel / mag if mag > 1e-9 else np.zeros(3)
    ut = v_rel_t / vt_mag if vt_mag > 1e-9 else np.zeros(3)
    return TriKinematics(n, c, triangle_area(tri), xr, v_rel, mag, up,
                         float(np.dot(up, n)), ut, alpha)


def viscous_drag(k: TriKinematics, cF, rho=RHO_WATER):
    fDrag = 0.5 * rho * cF * k.area * k.v_rel_mag
    vf = -k.ut * k.v_rel_mag
    f = vf * fDrag
    return f, np.cross(k.xr, f)


def pressure_drag(k: TriKinematics, params=DEFAULT_PARAMS):
    S = k.area
    v = k.v_rel_mag / params["vRDrag"]
    ct = k.cosTheta
    if ct >= 0.0:
        drag = -(params["cPDrag1"] * v + params["cPDrag2"] * v * v) * S * ct ** params["fPDrag"]
    else:
        drag = (params["cSDrag1"] * v + params["cSDrag2"] * v * v) * S * (-ct) ** params["fSDrag"]
    f = k.normal * drag
    return f, np.cross(k.xr, f)


def foil_lift(k: TriKinematics, AR, params=DEFAULT_PARAMS, rho=RHO_WATER,
              bottom_thresh=BOTTOM_THRESHOLD):
    """Port of ComputeFoilLiftForce *as implemented* (see audit finding F-03)."""
    nz = k.normal[2]
    if nz < bottom_thresh or k.v_rel_mag < 1e-4:
        return np.zeros(3), np.zeros(3)
    Cl_alpha = params["cLift1"] * 2.0 * math.pi
    a = k.alpha
    Cl = Cl_alpha * a if abs(a) < params["alphaStall"] else params["cLMax"] * (1.0 if a > 0 else -1.0)
    Cdi = Cl * Cl / (math.pi * AR + 1e-9)
    qA = 0.5 * rho * k.v_rel_mag ** 2 * k.area
    lift_dir = k.normal - k.up * float(np.dot(k.normal, k.up))
    ld = np.linalg.norm(lift_dir)
    if ld < 1e-9:
        return np.zeros(3), np.zeros(3)
    lift_dir /= ld
    F = lift_dir * (Cl * qA) + (-k.up) * (Cdi * qA)
    return F, np.cross(k.xr, F)


def aero_drag_face(normal, area_above, v_hull, v_wind, cd=1.0,
                   windward=True, rho=RHO_AIR):
    """Above-waterline face pressure force (Hydrodynamics.cc).

    ``windward=False`` reproduces the pre-audit face selection
    (vn = (v_wind - v_hull).n > 0, force along +n); ``windward=True`` the
    corrected one (vn = (v_hull - v_wind).n > 0, force along -n).  Both give
    the same net force on a body that is symmetric about the wind axis, but
    only the corrected one loads the faces that actually face the wind
    (Fossen Sec. 8.1, wind loads act on the projected windward area).
    """
    n = np.asarray(normal, dtype=float)
    if windward:
        vn = float(np.dot(np.asarray(v_hull) - np.asarray(v_wind), n))
        if vn <= 0.0:
            return np.zeros(3)
        return -(0.5 * rho * cd * area_above * vn * vn) * n
    vn = float(np.dot(np.asarray(v_wind) - np.asarray(v_hull), n))
    if vn <= 0.0:
        return np.zeros(3)
    return (0.5 * rho * cd * area_above * vn * vn) * n


# ---------------------------------------------------------------------------
# Simple closed hull meshes
# ---------------------------------------------------------------------------
def box_mesh(lx, ly, lz, center=(0.0, 0.0, 0.0)):
    """12-triangle closed box with outward normals, centred at ``center``."""
    cx, cy, cz = center
    x0, x1 = cx - lx / 2, cx + lx / 2
    y0, y1 = cy - ly / 2, cy + ly / 2
    z0, z1 = cz - lz / 2, cz + lz / 2
    P = np.array([[x0, y0, z0], [x1, y0, z0], [x1, y1, z0], [x0, y1, z0],
                  [x0, y0, z1], [x1, y0, z1], [x1, y1, z1], [x0, y1, z1]])
    faces = [
        (0, 2, 1), (0, 3, 2),   # bottom (n = -z)
        (4, 5, 6), (4, 6, 7),   # top    (n = +z)
        (0, 1, 5), (0, 5, 4),   # -y side
        (1, 2, 6), (1, 6, 5),   # +x side
        (2, 3, 7), (2, 7, 6),   # +y side
        (3, 0, 4), (3, 4, 7),   # -x side
    ]
    return [P[list(f)] for f in faces]


def transform_mesh(tris, R, t):
    return [np.asarray(tri) @ R.T + np.asarray(t) for tri in tris]


def signed_volume(tris):
    v = 0.0
    for a, b, c in tris:
        v += np.dot(a, np.cross(b, c)) / 6.0
    return v


def rotation_matrix(roll=0.0, pitch=0.0, yaw=0.0):
    """R_b^n(Theta) from Fossen Eq. 2.18 (zyx convention)."""
    cr, sr = math.cos(roll), math.sin(roll)
    cp, sp = math.cos(pitch), math.sin(pitch)
    cy, sy = math.cos(yaw), math.sin(yaw)
    return np.array([
        [cy * cp, -sy * cr + cy * sp * sr, sy * sr + cy * cr * sp],
        [sy * cp, cy * cr + sr * sp * sy, -cy * sr + sp * sy * cr],
        [-sp, cp * sr, cp * cr]])
