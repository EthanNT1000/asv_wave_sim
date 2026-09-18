"""Golden-value regression test for the geometry facade (CGAL removal).

Runs the ``geom_baseline`` tool built from ``gz-waves/test/geom`` (with
``-DBUILD_TESTING=ON``) and compares its JSON output against
``tests/geom/golden.json``. The golden file was captured in Phase 0 with the
CGAL backend and must hold, within tolerance, through every later phase.

Locate the tool with ``GEOM_BASELINE_BIN`` or pass ``--bin``; default is
``build/bin/geom_baseline`` relative to the repository root.

    python3 -m unittest tests/geom/test_geom_baseline.py -v
    python3 tests/geom/test_geom_baseline.py --update   # rewrite golden.json

Tolerances: relative 1e-9 / absolute 1e-7 on every number. A tolerance
change is not an acceptable way to make this test pass: if a phase moves a
value, the physics changed and the phase must be reported, not tuned.
"""
from __future__ import annotations

import json
import math
import os
import subprocess
import sys
import unittest

REPO = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", ".."))
GOLDEN = os.path.join(os.path.dirname(os.path.abspath(__file__)), "golden.json")
RTOL = 1e-9
ATOL = 1e-7


def baseline_bin():
    b = os.environ.get("GEOM_BASELINE_BIN")
    if b:
        return b
    return os.path.join(REPO, "build", "bin", "geom_baseline")


def run_baseline(extra=()):
    exe = baseline_bin()
    if not os.path.exists(exe):
        raise unittest.SkipTest(
            f"geom_baseline not found at {exe}; build gz-waves with "
            "-DBUILD_TESTING=ON and set GEOM_BASELINE_BIN")
    env = dict(os.environ)
    env.setdefault("GZ_PARTITION", "geom_baseline_test")  # isolate transport
    out = subprocess.run([exe, *extra], check=True, capture_output=True,
                         text=True, env=env).stdout
    return json.loads(out)


def compare(path, got, exp, errors):
    if isinstance(exp, dict):
        if not isinstance(got, dict):
            errors.append(f"{path}: expected object"); return
        for k in exp:
            if k not in got:
                errors.append(f"{path}.{k}: missing"); continue
            compare(f"{path}.{k}", got[k], exp[k], errors)
        for k in got:
            if k not in exp and k != "bench":
                errors.append(f"{path}.{k}: unexpected key")
    elif isinstance(exp, list):
        if not isinstance(got, list) or len(got) != len(exp):
            errors.append(f"{path}: length {len(got) if isinstance(got, list) else '?'} != {len(exp)}")
            return
        for i, (g, e) in enumerate(zip(got, exp)):
            compare(f"{path}[{i}]", g, e, errors)
    elif exp is None:
        if got is not None:
            errors.append(f"{path}: expected null, got {got}")
    elif isinstance(exp, bool):
        if got != exp:
            errors.append(f"{path}: {got} != {exp}")
    else:
        if got is None or not isinstance(got, (int, float)):
            errors.append(f"{path}: expected number, got {got!r}"); return
        if not math.isclose(got, exp, rel_tol=RTOL, abs_tol=ATOL):
            errors.append(f"{path}: {got!r} != {exp!r} (diff {got - exp:.3e})")


class GeomBaselineTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.got = run_baseline()
        with open(GOLDEN) as f:
            cls.exp = json.load(f)

    def _check(self, key):
        errors = []
        compare(key, self.got[key], self.exp[key], errors)
        self.assertFalse(errors, "\n".join(errors[:40]))

    def test_wave_height_at_sample_points(self):
        self._check("wave_height")

    def test_hydrodynamics_force_torque_volume(self):
        self._check("hydrodynamics")

    def test_ray_mesh_first_intersection(self):
        self._check("ray_mesh")

    def test_still_water_force_is_buoyancy(self):
        """Sanity on the golden data itself: a 10x4x2 box half submerged."""
        rest = self.got["hydrodynamics"]["still_box_10x4x2_rest"]
        self.assertAlmostEqual(rest["displaced_volume_from_fz"], 10 * 4 * 1.0, places=6)
        self.assertAlmostEqual(rest["force"][0], 0.0, places=6)
        self.assertAlmostEqual(rest["force"][1], 0.0, places=6)


if __name__ == "__main__":
    if "--update" in sys.argv:
        data = run_baseline()
        data.pop("bench", None)
        with open(GOLDEN, "w") as f:
            json.dump(data, f, indent=1)
        print(f"wrote {GOLDEN}")
    else:
        unittest.main(verbosity=2)
