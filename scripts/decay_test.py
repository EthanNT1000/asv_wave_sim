#!/usr/bin/env python3
"""Free-decay test (roll, pitch, heave) for the gz-waves hydrodynamics model.

Integrates a rigid hull with the numpy port of ``Physics.cc``
(``tests/hydro/hydro_reference.py``): exact per-triangle hydrostatics,
ITTC-57 viscous drag, Kerner pressure drag, foil lift (as implemented) and
the per-DOF Fossen damping, with the same fixed-step explicit force
evaluation as the Gazebo plugin.  For each DOF it releases the hull from an
initial offset in still water and reports the natural period and damping
ratio (log decrement), plus the linear estimate ``2 pi sqrt(M_ii / C_ii)``
from the hydrostatic stiffness.

Examples
--------
    python3 scripts/decay_test.py                       # WAM-V-like twin pontoon
    python3 scripts/decay_test.py --hull box            # 4 x 1.5 x 0.6 m, 1000 kg
    python3 scripts/decay_test.py --com-offset 0.3 0 0.1 --legacy-kinematics
    python3 scripts/decay_test.py --added-mass 10 90 180 40 150 80
    python3 scripts/decay_test.py --json out.json

``--legacy-kinematics`` reproduces the pre-audit behaviour in which the link
origin velocity was used as the CoM velocity (audit finding F-05); it only
differs from the fixed behaviour when ``--com-offset`` is non-zero.
``--added-mass`` adds a diagonal M_A to the mass matrix used by the
integrator, emulating an SDF ``<fluid_added_mass>`` (audit finding F-01).
"""
from __future__ import annotations

import argparse
import json
import math
import os
import sys

import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)),
                                "..", "tests", "hydro"))
import hydro_reference as hr  # noqa: E402


# ---------------------------------------------------------------------------
# Hulls
# ---------------------------------------------------------------------------
def make_hull(name):
    """Return (triangles in body frame with origin at the keel plane centre,
    mass, inertia_cg, hull length, waterplane area, KB-independent data)."""
    if name == "twin":
        # WAM-V-like: two box pontoons, inertia from gz-waves-models/wam-v.
        L, b, h, sep = 4.9, 0.30, 0.50, 2.44
        tris = (hr.box_mesh(L, b, h, center=(0.0, +sep / 2, h / 2))
                + hr.box_mesh(L, b, h, center=(0.0, -sep / 2, h / 2)))
        m = 180.0
        I = np.diag([120.0, 393.0, 446.0])
        Awp = 2 * L * b
        I_wp_T = 2 * (L * b ** 3 / 12 + L * b * (sep / 2) ** 2)
        I_wp_L = 2 * b * L ** 3 / 12
        return tris, m, I, L, Awp, I_wp_T, I_wp_L, h
    if name == "box":
        L, B, H = 4.0, 1.5, 0.6
        m = 1000.0
        I = m / 12 * np.diag([B ** 2 + H ** 2, L ** 2 + H ** 2, L ** 2 + B ** 2])
        return (hr.box_mesh(L, B, H, center=(0, 0, H / 2)), m, I, L, L * B,
                L * B ** 3 / 12, B * L ** 3 / 12, H)
    raise ValueError(name)


# ---------------------------------------------------------------------------
# Rigid-body + hydrodynamics integrator (mirrors the plugin's force pipeline)
# ---------------------------------------------------------------------------
class Sim:
    def __init__(self, hull, com_body, params, dt=1e-3, added_mass=None,
                 foil=True, legacy_kinematics=False, r_g=(0.0, 0.0, 0.0)):
        (self.tris_body, self.m, self.I_cg, self.L, self.Awp,
         self.Iwp_T, self.Iwp_L, self.h) = make_hull(hull)
        # Mesh vertices relative to the CoM (body frame); com_body is the CoM
        # position in the mesh (keel-plane) frame, r_g the CoM position
        # relative to the SDF link origin (the inertial <pose>).
        self.com_body = np.asarray(com_body, dtype=float)
        self.r_g = np.asarray(r_g, dtype=float)
        self.drag = True
        self.tris_cm = [t - self.com_body for t in self.tris_body]
        self.params = params
        self.dt = dt
        self.foil = foil
        self.legacy = legacy_kinematics
        self.M_A = np.zeros(6) if added_mass is None else np.asarray(added_mass, float)
        self.m_eff = self.m + self.M_A[:3]                 # per-axis (diagonal M_A)
        self.I_eff_b = self.I_cg + np.diag(self.M_A[3:])
        self.depth = lambda p: -p[2]

    # hydrostatic equilibrium draft for the upright hull (bisection)
    def equilibrium_z(self):
        lo, hi = -self.h, self.h
        for _ in range(80):
            mid = 0.5 * (lo + hi)
            F, _, _ = hr.hydrostatic_wrench(
                hr.transform_mesh(self.tris_cm, np.eye(3), [0, 0, mid]),
                self.depth, np.array([0, 0, mid]))
            if F[2] + self.m * hr.GRAVITY > 0:      # too much buoyancy -> raise
                lo = mid
            else:
                hi = mid
        return 0.5 * (lo + hi)

    def wrench(self, x, R, v_com, omega):
        """Hydrodynamic force/torque about the CoM, vectorised over the
        submerged sub-triangles (same terms as Hydrodynamics::Update)."""
        tris = np.asarray(hr.transform_mesh(self.tris_cm, R, x))
        F, tau, sub = hr.hydrostatic_wrench_fast(tris, x)
        # Pre-audit behaviour: the plugin used the link-origin velocity as the
        # CoM velocity.  v_link = v_com - omega x (R r_g).
        v_used = v_com - np.cross(omega, R @ self.r_g) if self.legacy else v_com
        if len(sub) == 0 or not self.drag:
            return F, tau
        p = self.params
        rho = hr.RHO_WATER
        e1 = sub[:, 1] - sub[:, 0]
        e2 = sub[:, 2] - sub[:, 0]
        nA2 = np.cross(e1, e2)
        A = 0.5 * np.linalg.norm(nA2, axis=1)
        n = nA2 / np.maximum(2 * A, 1e-300)[:, None]
        c = sub.mean(axis=1)
        xr = c - x
        v_rel = v_used + np.cross(omega, xr)              # still water: v_fluid = 0
        mag = np.linalg.norm(v_rel, axis=1)
        vdn = np.einsum("ij,ij->i", v_rel, n)
        v_t = v_rel - n * vdn[:, None]
        vt_mag = np.linalg.norm(v_t, axis=1)
        up = np.where(mag[:, None] > 1e-9, v_rel / np.maximum(mag, 1e-300)[:, None], 0.0)
        ut = np.where(vt_mag[:, None] > 1e-9, v_t / np.maximum(vt_mag, 1e-300)[:, None], 0.0)
        cosT = np.einsum("ij,ij->i", up, n)

        # viscous drag (ITTC-57, Reynolds number from the CoM speed and L)
        Rn = np.linalg.norm(v_used) * self.L / hr.NU_WATER
        cF = hr.ittc57_cf(Rn)
        f_visc = -ut * (0.5 * rho * cF * A * mag * mag)[:, None]

        # pressure drag
        v = mag / p["vRDrag"]
        pos = cosT >= 0.0
        drag = np.where(
            pos,
            -(p["cPDrag1"] * v + p["cPDrag2"] * v * v) * A * np.abs(cosT) ** p["fPDrag"],
            (p["cSDrag1"] * v + p["cSDrag2"] * v * v) * A * np.abs(cosT) ** p["fSDrag"])
        f_pres = n * drag[:, None]

        f_all = f_visc + f_pres

        # foil lift as implemented (gate n_z > BOTTOM_THRESHOLD, force along +n)
        if self.foil:
            bottom = n[:, 2] > hr.BOTTOM_THRESHOLD
            span = max(np.ptp(sub[:, :, 1]), 1e-9)
            AR = max(span * span / (A[bottom].sum() + 1e-9), 0.5)
            alpha = np.arctan2(vdn, vt_mag + 1e-9)
            Cl_alpha = p["cLift1"] * 2.0 * math.pi
            Cl = np.where(np.abs(alpha) < p["alphaStall"], Cl_alpha * alpha,
                          p["cLMax"] * np.sign(alpha))
            Cdi = Cl * Cl / (math.pi * AR + 1e-9)
            qA = 0.5 * rho * mag * mag * A
            lift_dir = n - up * np.einsum("ij,ij->i", n, up)[:, None]
            ld = np.linalg.norm(lift_dir, axis=1)
            ok = bottom & (mag >= 1e-4) & (ld >= 1e-9)
            lift_dir = np.where(ok[:, None], lift_dir / np.maximum(ld, 1e-300)[:, None], 0.0)
            f_foil = lift_dir * (Cl * qA)[:, None] - up * (Cdi * qA)[:, None]
            f_all = f_all + np.where(ok[:, None], f_foil, 0.0)

        F = F + f_all.sum(axis=0)
        tau = tau + np.cross(xr, f_all).sum(axis=0)
        f, tq = hr.damping_force_world(R, v_used, omega, np.zeros(3), p)
        return F + f, tau + tq

    def run(self, x0, R0, t_end):
        x, R = np.array(x0, float), np.array(R0, float)
        v, w = np.zeros(3), np.zeros(3)
        n = int(round(t_end / self.dt))
        out = np.zeros((n, 4))
        g = np.array([0, 0, hr.GRAVITY])
        for i in range(n):
            F, tau = self.wrench(x, R, v, w)
            # forces are explicit (previous state), rigid body semi-implicit
            a = (F + self.m * g) / self.m_eff
            I_w = R @ self.I_eff_b @ R.T
            alpha = np.linalg.solve(I_w, tau - np.cross(w, I_w @ w))
            v = v + self.dt * a
            w = w + self.dt * alpha
            x = x + self.dt * v
            R = _expm_so3(w * self.dt) @ R
            roll = math.atan2(R[2, 1], R[2, 2])
            pitch = -math.asin(max(-1.0, min(1.0, R[2, 0])))
            out[i] = (i * self.dt, x[2], roll, pitch)
        return out


def _expm_so3(phi):
    th = np.linalg.norm(phi)
    if th < 1e-12:
        return np.eye(3)
    K = hr.skew(phi / th)
    return np.eye(3) + math.sin(th) * K + (1 - math.cos(th)) * K @ K


# ---------------------------------------------------------------------------
# Decay analysis
# ---------------------------------------------------------------------------
def analyse(t, s):
    """Natural period from up-crossings, damping ratio from log decrement."""
    s = s - s[-len(s) // 5:].mean()            # remove residual offset
    # analyse only while the motion is above 2 % of its initial amplitude,
    # otherwise numerical noise in the settled tail biases the crossings
    a0 = np.abs(s[:max(1, len(s) // 50)]).max()
    alive = np.where(np.abs(s) > 0.02 * a0)[0]
    n_end = int(alive[-1]) + 1 if len(alive) else len(s)
    t, s = t[:n_end], s[:n_end]
    up = np.where((s[:-1] < 0) & (s[1:] >= 0))[0]
    peaks = [i for i in range(1, len(s) - 1) if s[i] > s[i - 1] and s[i] >= s[i + 1] and s[i] > 0]
    T = float(np.mean(np.diff(t[up]))) if len(up) >= 2 else float("nan")
    zeta = float("nan")
    if len(peaks) >= 3:
        A = s[peaks[:6]]
        d = np.log(A[:-1] / A[1:])
        d = d[np.isfinite(d) & (d > 0)]
        if len(d):
            delta = float(np.mean(d))
            zeta = delta / math.sqrt(4 * math.pi ** 2 + delta ** 2)
    return T, zeta, len(peaks)


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--hull", choices=["twin", "box"], default="twin")
    ap.add_argument("--dt", type=float, default=1e-3, help="step size (Gazebo default 1 ms)")
    ap.add_argument("--t-end", type=float, default=12.0)
    ap.add_argument("--com-offset", type=float, nargs=3, default=(0.0, 0.0, 0.0),
                    metavar=("DX", "DY", "DZ"),
                    help="CoM offset from the hull reference (link) origin, body frame [m]")
    ap.add_argument("--kg", type=float, default=None,
                    help="CoM height above the keel [m] (default: 0.6 h twin, 0.5 H box)")
    ap.add_argument("--added-mass", type=float, nargs=6, default=None,
                    metavar=("A11", "A22", "A33", "A44", "A55", "A66"))
    ap.add_argument("--legacy-kinematics", action="store_true")
    ap.add_argument("--no-foil", action="store_true")
    ap.add_argument("--no-drag", action="store_true",
                    help="disable viscous/pressure drag and per-DOF damping (hydrostatics only)")
    ap.add_argument("--heave0", type=float, default=0.03, help="initial heave offset [m]")
    ap.add_argument("--roll0", type=float, default=math.radians(2.0),
                    help="initial roll [rad]; keep both pontoons wet for the twin hull")
    ap.add_argument("--pitch0", type=float, default=math.radians(1.0), help="initial pitch [rad]")
    ap.add_argument("--json", default=None)
    args = ap.parse_args()

    params = dict(hr.DEFAULT_PARAMS)
    _, _, _, _, _, _, _, h = make_hull(args.hull)
    kg = args.kg if args.kg is not None else (0.6 * h if args.hull == "twin" else 0.5 * h)
    com_body = np.array(args.com_offset, float) + np.array([0.0, 0.0, kg])

    sim = Sim(args.hull, com_body, params, dt=args.dt, added_mass=args.added_mass,
              foil=not args.no_foil, legacy_kinematics=args.legacy_kinematics,
              r_g=args.com_offset)
    sim.drag = not args.no_drag
    z_eq = sim.equilibrium_z()
    kg = float(com_body[2])      # CoM height above the keel plane
    T_draft = kg - z_eq          # keel below the surface (CoM at z_eq, keel at z_eq - kg)

    # linear estimates (Fossen Sec. 4.2: C33 = rho g Awp, C44 = rho g V GM_T, C55 = rho g V GM_L)
    rho_g = hr.RHO_WATER * 9.81
    V = sim.m / hr.RHO_WATER
    KB = T_draft / 2
    GM_T = KB + sim.Iwp_T / V - kg          # box pontoons: exact for small angles
    GM_L = KB + sim.Iwp_L / V - kg
    C = {"heave": rho_g * sim.Awp, "roll": rho_g * V * GM_T, "pitch": rho_g * V * GM_L}
    Mii = {"heave": sim.m_eff[2], "roll": sim.I_eff_b[0, 0], "pitch": sim.I_eff_b[1, 1]}

    results = {}
    cases = {
        "heave": (np.array([0, 0, z_eq + args.heave0]), np.eye(3), 1),
        "roll": (np.array([0, 0, z_eq]), hr.rotation_matrix(roll=args.roll0), 2),
        "pitch": (np.array([0, 0, z_eq]), hr.rotation_matrix(pitch=args.pitch0), 3),
    }
    print(f"hull={args.hull} m={sim.m:.0f} kg draft={T_draft:.3f} m KG={kg:.3f} m "
          f"GM_T={GM_T:.3f} m GM_L={GM_L:.3f} m dt={args.dt} "
          f"com_offset={tuple(args.com_offset)} legacy={args.legacy_kinematics} "
          f"foil={not args.no_foil} M_A={args.added_mass}")
    print(f"{'DOF':6s} {'T_n [s]':>9s} {'zeta':>8s} {'T_lin [s]':>10s} {'peaks':>6s}")
    for dof, (x0, R0, col) in cases.items():
        out = sim.run(x0, R0, args.t_end)
        T, zeta, npk = analyse(out[:, 0], out[:, col])
        T_lin = 2 * math.pi * math.sqrt(Mii[dof] / C[dof]) if C[dof] > 0 else float("nan")
        results[dof] = {"T_n": T, "zeta": zeta, "T_linear": T_lin, "peaks": npk}
        print(f"{dof:6s} {T:9.3f} {zeta:8.4f} {T_lin:10.3f} {npk:6d}")
    if args.json:
        with open(args.json, "w") as f:
            json.dump({"args": vars(args), "draft": T_draft, "GM_T": GM_T, "GM_L": GM_L,
                       "results": results}, f, indent=2)
    return 0


if __name__ == "__main__":
    sys.exit(main())
