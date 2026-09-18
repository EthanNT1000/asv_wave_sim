"""Structural property tests for the gz-waves hydrodynamics model.

Run with either

    python3 -m unittest tests/hydro/test_properties.py -v
    python3 -m pytest tests/hydro -v

The tests use the pure-numpy port in ``hydro_reference.py`` (mirrors
``gz-waves/src/Physics.cc``) plus the Fossen matrix constructions, and read
the WAM-V model SDF for the rigid-body inertia. They need only numpy.

Properties asserted (see docs/hydrodynamics_audit.md):
  * M = M_RB + M_A symmetric positive definite
  * C(nu) derived from M is skew-symmetric (Fossen Theorem 3.2), incl. C_A
  * D(nu) = D_L + D_n(nu) positive definite, dissipative
  * zero net force and moment at rest in still water (buoyancy = weight)
  * every relative-velocity term is dissipative (F . v_rel <= 0)
  * restoring stiffness matches rho g V GM_T / GM_L (Fossen Ch. 4)
"""
from __future__ import annotations

import math
import os
import sys
import unittest
import xml.etree.ElementTree as ET

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import hydro_reference as hr  # noqa: E402

REPO = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", ".."))
WAMV_SDF = os.path.join(REPO, "gz-waves-models", "models", "wam-v", "model.sdf")


def read_sdf_inertial(path=WAMV_SDF, link="base_link"):
    """Return (m, I_cg 3x3, r_g, M_A 6x6 or None) from an SDF <inertial>."""
    root = ET.parse(path).getroot()
    for lk in root.iter("link"):
        if lk.get("name") != link:
            continue
        inertial = lk.find("inertial")
        m = float(inertial.findtext("mass"))
        I = inertial.find("inertia")
        ixx, iyy, izz = (float(I.findtext(k)) for k in ("ixx", "iyy", "izz"))
        ixy, ixz, iyz = (float(I.findtext(k) or 0.0) for k in ("ixy", "ixz", "iyz"))
        I_cg = np.array([[ixx, ixy, ixz], [ixy, iyy, iyz], [ixz, iyz, izz]])
        pose = inertial.findtext("pose")
        r_g = np.array([float(v) for v in pose.split()[:3]]) if pose else np.zeros(3)
        M_A = None
        fam = inertial.find("fluid_added_mass")
        if fam is not None:
            keys = ["xx", "xy", "xz", "xp", "xq", "xr", "yy", "yz", "yp", "yq", "yr",
                    "zz", "zp", "zq", "zr", "pp", "pq", "pr", "qq", "qr", "rr"]
            vals = {k: float(fam.findtext(k) or 0.0) for k in keys}
            names = "xyzpqr"
            M_A = np.zeros((6, 6))
            for i in range(6):
                for j in range(i, 6):
                    k = names[i] + names[j]
                    M_A[i, j] = M_A[j, i] = vals.get(k, 0.0)
        return m, I_cg, r_g, M_A
    raise KeyError(link)


def read_sdf_hydro_params(path=WAMV_SDF):
    params = dict(hr.DEFAULT_PARAMS)
    root = ET.parse(path).getroot()
    hydro = next(root.iter("hydrodynamics"), None)
    if hydro is not None:
        for child in hydro:
            if child.tag in params:
                params[child.tag] = float(child.text)
    return params


# Representative added mass for a WAM-V-class twin-pontoon hull, used only to
# exercise the C_A(nu_r) derivation (the SDF carries no <fluid_added_mass>,
# audit finding F-01).
M_A_EXAMPLE = np.diag([10.0, 90.0, 180.0, 40.0, 150.0, 80.0])


class MassMatrixTests(unittest.TestCase):
    def setUp(self):
        self.m, self.I_cg, self.r_g, self.M_A_sdf = read_sdf_inertial()
        self.M_RB = hr.rigid_body_mass_matrix(self.m, self.I_cg, self.r_g)
        self.M_A = self.M_A_sdf if self.M_A_sdf is not None else M_A_EXAMPLE
        self.M = self.M_RB + self.M_A

    def test_mass_matrix_spd(self):
        for name, M in (("M_RB", self.M_RB), ("M_A", self.M_A), ("M", self.M)):
            with self.subTest(matrix=name):
                np.testing.assert_allclose(M, M.T, atol=1e-12, err_msg=f"{name} not symmetric")
                eig = np.linalg.eigvalsh(M)
                self.assertGreater(eig.min(), 0.0, f"{name} not positive definite: {eig}")

    def test_rigid_body_mass_matrix_with_cg_offset_is_spd(self):
        M = hr.rigid_body_mass_matrix(self.m, self.I_cg, r_g=(0.3, -0.1, 0.2))
        np.testing.assert_allclose(M, M.T, atol=1e-12)
        self.assertGreater(np.linalg.eigvalsh(M).min(), 0.0)


class CoriolisTests(unittest.TestCase):
    """C(nu) built from M by Fossen Theorem 3.2 must be skew-symmetric."""

    def setUp(self):
        self.rng = np.random.default_rng(42)
        m, I_cg, r_g, M_A_sdf = read_sdf_inertial()
        self.M_RB = hr.rigid_body_mass_matrix(m, I_cg, r_g)
        self.M_A = M_A_sdf if M_A_sdf is not None else M_A_EXAMPLE

    def test_C_RB_skew_symmetric(self):
        for _ in range(50):
            nu = self.rng.normal(size=6)
            C = hr.coriolis_from_mass_matrix(self.M_RB, nu)
            np.testing.assert_allclose(C + C.T, 0.0, atol=1e-9)
            self.assertAlmostEqual(float(nu @ C @ nu), 0.0, places=8)   # no work

    def test_C_A_skew_symmetric(self):
        for _ in range(50):
            nu_r = self.rng.normal(size=6)
            C = hr.coriolis_from_mass_matrix(self.M_A, nu_r)
            np.testing.assert_allclose(C + C.T, 0.0, atol=1e-9)

    def test_C_RB_with_cg_offset_skew_symmetric(self):
        m, I_cg, _, _ = read_sdf_inertial()
        M = hr.rigid_body_mass_matrix(m, I_cg, r_g=(0.3, -0.1, 0.2))
        for _ in range(20):
            nu = self.rng.normal(size=6)
            C = hr.coriolis_from_mass_matrix(M, nu)
            np.testing.assert_allclose(C + C.T, 0.0, atol=1e-9)

    def test_munk_moment_follows_from_M_A(self):
        """Fossen Sec. 6.3: for u_r > 0, v_r > 0 the yaw Munk moment is
        (Y_vdot - X_udot) u_r v_r, i.e. -(A22 - A11) u v in the derived C_A."""
        u, v = 2.0, 0.5
        nu_r = np.array([u, v, 0.0, 0.0, 0.0, 0.0])
        N = hr.munk_moment(self.M_A, nu_r)[2]
        X_udot, Y_vdot = -self.M_A[0, 0], -self.M_A[1, 1]
        self.assertAlmostEqual(N, (Y_vdot - X_udot) * u * v, places=9)
        # |Y_vdot| > |X_udot|: bow turns away from the inflow (destabilising)
        self.assertLess(N, 0.0)


class DampingTests(unittest.TestCase):
    def setUp(self):
        self.rng = np.random.default_rng(7)
        self.params = read_sdf_hydro_params()

    def test_damping_positive_definite(self):
        for _ in range(50):
            nu = self.rng.normal(size=6) * 3.0
            D = hr.damping_matrix(nu, self.params)
            np.testing.assert_allclose(D, D.T, atol=0.0)
            self.assertGreater(np.linalg.eigvalsh(D).min(), 0.0)
            self.assertGreaterEqual(float(nu @ D @ nu), 0.0)          # dissipative

    def test_damping_is_D_L_plus_D_n(self):
        nu = np.array([1.0, -2.0, 0.5, 0.1, -0.2, 0.3])
        D = hr.damping_matrix(nu, self.params)
        D_L = hr.damping_matrix(np.zeros(6), self.params)
        D_n = D - D_L
        self.assertTrue(np.all(np.diag(D_n) >= 0.0))
        for i, key in enumerate("UVWPQN"):
            self.assertAlmostEqual(D_n[i, i], self.params[f"cDamp{key}2"] * abs(nu[i]))

    def test_damping_world_frame_dissipative(self):
        for _ in range(20):
            R = hr.rotation_matrix(*self.rng.uniform(-1, 1, size=3))
            v, w, vc = (self.rng.normal(size=3) for _ in range(3))
            vc[2] = 0.0
            f, tau = hr.damping_force_world(R, v, w, vc, self.params)
            power = float(np.dot(f, v - vc) + np.dot(tau, w))
            self.assertLessEqual(power, 1e-12)


class StillWaterTests(unittest.TestCase):
    """Zero net force/moment at rest in still water and exact hydrostatics."""

    L, B, H, m = 4.0, 1.5, 0.6, 1000.0

    def setUp(self):
        self.T = self.m / (hr.RHO_WATER * self.L * self.B)      # equilibrium draft
        self.com = np.array([0.0, 0.0, self.H / 2 - self.T])   # KG = H/2
        self.depth = lambda p: -p[2]                             # flat surface z = 0
        self.tris = hr.transform_mesh(hr.box_mesh(self.L, self.B, self.H),
                                      np.eye(3), self.com)

    def test_zero_net_force_at_rest(self):
        F, tau, sub = hr.hydrostatic_wrench(self.tris, self.depth, self.com)
        weight = np.array([0.0, 0.0, self.m * hr.GRAVITY])
        np.testing.assert_allclose(F + weight, 0.0, atol=1e-6 * self.m * 9.81)
        np.testing.assert_allclose(tau, 0.0, atol=1e-6)
        self.assertGreater(len(sub), 0)
        # velocity-dependent terms vanish identically at nu = 0
        R = np.eye(3)
        f, t = hr.damping_force_world(R, np.zeros(3), np.zeros(3), np.zeros(3))
        np.testing.assert_allclose(f, 0.0)
        np.testing.assert_allclose(t, 0.0)
        for tri in sub:
            k = hr.point_kinematics(tri, self.com, np.zeros(3), np.zeros(3),
                                    lambda c: np.zeros(3))
            self.assertEqual(k.v_rel_mag, 0.0)
            np.testing.assert_allclose(hr.viscous_drag(k, 0.003)[0], 0.0)
            np.testing.assert_allclose(hr.pressure_drag(k)[0], 0.0)
            np.testing.assert_allclose(hr.foil_lift(k, AR=2.0)[0], 0.0)

    def test_submerged_body_force_is_rho_g_V(self):
        c = np.array([0.0, 0.0, -3.0])
        R = hr.rotation_matrix(0.3, 0.5, 0.2)
        tris = hr.transform_mesh(hr.box_mesh(1.0, 1.0, 1.0), R, c)
        F, tau, _ = hr.hydrostatic_wrench(tris, self.depth, c)
        np.testing.assert_allclose(F, [0.0, 0.0, hr.RHO_WATER * 9.81], rtol=1e-9, atol=1e-6)
        np.testing.assert_allclose(tau, 0.0, atol=1e-6)      # CoP formulas exact

    def test_vectorised_hydrostatics_matches_cop_port(self):
        """hydrostatic_wrench_fast (exact triangle rule) must equal the
        centre-of-pressure construction ported from Physics.cc."""
        rng = np.random.default_rng(11)
        for _ in range(10):
            R = hr.rotation_matrix(*rng.uniform(-0.6, 0.6, size=3))
            c = rng.normal(size=3) * 0.3
            tris = hr.transform_mesh(hr.box_mesh(self.L, self.B, self.H), R, c)
            F1, t1, _ = hr.hydrostatic_wrench(tris, self.depth, c)
            F2, t2, _ = hr.hydrostatic_wrench_fast(tris, c)
            np.testing.assert_allclose(F1, F2, rtol=1e-10, atol=1e-8)
            np.testing.assert_allclose(t1, t2, rtol=1e-10, atol=1e-8)

    def test_restoring_moments_match_GM(self):
        """g(eta) linearises to rho g V GM_T phi, rho g V GM_L theta (Fossen Sec. 4.2)."""
        V = self.L * self.B * self.T
        KB, KG = self.T / 2, self.H / 2
        GM_T = KB + (self.L * self.B ** 3 / 12) / V - KG
        GM_L = KB + (self.B * self.L ** 3 / 12) / V - KG
        phi, theta = 0.02, 0.01
        tris = hr.transform_mesh(hr.box_mesh(self.L, self.B, self.H),
                                 hr.rotation_matrix(roll=phi), self.com)
        _, tau, _ = hr.hydrostatic_wrench(tris, self.depth, self.com)
        self.assertAlmostEqual(tau[0], -hr.RHO_WATER * 9.81 * V * GM_T * phi, delta=1.0)
        tris = hr.transform_mesh(hr.box_mesh(self.L, self.B, self.H),
                                 hr.rotation_matrix(pitch=theta), self.com)
        _, tau, _ = hr.hydrostatic_wrench(tris, self.depth, self.com)
        self.assertAlmostEqual(tau[1], -hr.RHO_WATER * 9.81 * V * GM_L * theta, delta=1.0)
        self.assertGreater(GM_T, 0.0)
        self.assertGreater(GM_L, 0.0)


class DissipationTests(unittest.TestCase):
    """Every relative-velocity force term must do non-positive work on v_rel."""

    def setUp(self):
        self.rng = np.random.default_rng(3)

    def _random_tri(self):
        return self.rng.normal(size=(3, 3))

    def test_viscous_and_pressure_drag_dissipative(self):
        for _ in range(200):
            tri = self._random_tri()
            com = self.rng.normal(size=3)
            v, w = self.rng.normal(size=3), self.rng.normal(size=3)
            vf = self.rng.normal(size=3)
            k = hr.point_kinematics(tri, com, v, w, lambda c: vf)
            fv, _ = hr.viscous_drag(k, hr.ittc57_cf(1e6))
            fp, _ = hr.pressure_drag(k)
            self.assertLessEqual(float(np.dot(fv, k.v_rel)), 1e-9)
            self.assertLessEqual(float(np.dot(fp, k.v_rel)), 1e-9)

    def test_foil_induced_drag_dissipative_lift_does_no_work(self):
        for _ in range(200):
            tri = self._random_tri()
            k = hr.point_kinematics(tri, np.zeros(3), self.rng.normal(size=3),
                                    self.rng.normal(size=3), lambda c: np.zeros(3))
            F, _ = hr.foil_lift(k, AR=2.0, bottom_thresh=-10.0)
            self.assertLessEqual(float(np.dot(F, k.v_rel)), 1e-9)

    def test_wind_load_dissipative_and_on_windward_faces(self):
        for _ in range(100):
            n = self.rng.normal(size=3)
            n /= np.linalg.norm(n)
            v_hull, v_wind = self.rng.normal(size=3), self.rng.normal(size=3)
            f = hr.aero_drag_face(n, 1.0, v_hull, v_wind, windward=True)
            self.assertLessEqual(float(np.dot(f, v_hull - v_wind)), 1e-9)
            if np.linalg.norm(f) > 0:
                self.assertGreater(float(np.dot(v_hull - v_wind, n)), 0.0)  # windward
        # head wind on a box: force on the +x (bow) face, not on the -x face
        bow = hr.aero_drag_face([1, 0, 0], 1.0, [0, 0, 0], [-10, 0, 0], windward=True)
        stern = hr.aero_drag_face([-1, 0, 0], 1.0, [0, 0, 0], [-10, 0, 0], windward=True)
        np.testing.assert_allclose(bow, [-0.5 * hr.RHO_AIR * 100, 0, 0])
        np.testing.assert_allclose(stern, 0.0)


class FoilLiftReferenceTests(unittest.TestCase):
    """Audit finding F-02 (not applied): lift on a bow-up planing bottom must
    be upward and the bottom gate must accept downward-facing outward normals.
    Expected to fail against the current formula; remove the decorator once the
    fix in docs/hydrodynamics_audit.md is applied."""

    @unittest.expectedFailure
    def test_foil_lift_sign_reference(self):
        tau_trim = 0.05
        n = np.array([math.sin(tau_trim), 0.0, -math.cos(tau_trim)])   # outward, bottom
        u = np.cross(n, [0.0, 1.0, 0.0]); u /= np.linalg.norm(u)
        v = np.cross(n, u)
        tri = np.array([np.zeros(3), u, v])
        if np.dot(hr.triangle_normal(tri), n) < 0:
            tri = tri[[0, 2, 1]]
        k = hr.point_kinematics(tri, np.zeros(3), np.array([3.0, 0.0, 0.0]),
                                np.zeros(3), lambda c: np.zeros(3))
        F, _ = hr.foil_lift(k, AR=2.0)      # current gate: nz > 0.3 rejects the bottom
        self.assertGreater(np.linalg.norm(F), 0.0, "bottom face rejected by gate")
        self.assertGreater(F[2], 0.0, "planing lift must be upward")


if __name__ == "__main__":
    unittest.main(verbosity=2)
