"""Golden-value regression test for the FFT backend (FFTW -> PocketFFT).

Runs the ``fft_baseline`` tool built from ``gz-waves/test/fft`` (with
``-DBUILD_TESTING=ON``) and compares its JSON output against
``tests/fft/golden.json``. The golden file was captured with the FFTW3
backend; the PocketFFT backend must reproduce every value within the
tolerance below.

Locate the tool with ``FFT_BASELINE_BIN`` or pass ``--bin``; default is
``build/bin/fft_baseline`` relative to the repository root.

    python3 -m unittest tests/fft/test_fft_baseline.py -v
    python3 tests/fft/test_fft_baseline.py --update   # rewrite golden.json

Tolerance: absolute 1e-10 on every number (all transforms are double
precision; the single-precision tolerance of 1e-5 is not exercised because
the library has no float transforms). A tolerance change is not an
acceptable way to make this test pass.
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
ATOL = 1e-10


def baseline_bin():
    b = os.environ.get("FFT_BASELINE_BIN")
    if b:
        return b
    return os.path.join(REPO, "build", "bin", "fft_baseline")


def run_baseline(extra=()):
    exe = baseline_bin()
    if not os.path.exists(exe):
        raise unittest.SkipTest(
            f"fft_baseline not found at {exe}; build gz-waves with "
            "-DBUILD_TESTING=ON and set FFT_BASELINE_BIN")
    out = subprocess.run([exe, *extra], check=True, capture_output=True,
                         text=True).stdout
    return json.loads(out)


def compare(path, got, exp, errors, stats):
    if isinstance(exp, dict):
        if not isinstance(got, dict):
            errors.append(f"{path}: expected object"); return
        for k in exp:
            if k not in got:
                errors.append(f"{path}.{k}: missing"); continue
            compare(f"{path}.{k}", got[k], exp[k], errors, stats)
        for k in got:
            if k not in exp and k != "bench":
                errors.append(f"{path}.{k}: unexpected key")
    elif isinstance(exp, list):
        if not isinstance(got, list) or len(got) != len(exp):
            errors.append(f"{path}: length mismatch"); return
        for i, (g, e) in enumerate(zip(got, exp)):
            compare(f"{path}[{i}]", g, e, errors, stats)
    elif isinstance(exp, str):
        pass  # backend name is informational
    elif exp is None:
        if got is not None:
            errors.append(f"{path}: expected null, got {got}")
    else:
        if got is None or not isinstance(got, (int, float)):
            errors.append(f"{path}: expected number, got {got!r}"); return
        d = abs(got - exp)
        stats["n"] += 1
        stats["max_abs_diff"] = max(stats["max_abs_diff"], d)
        if d > ATOL:
            errors.append(f"{path}: {got!r} != {exp!r} (diff {got - exp:.3e})")


class FftBaselineTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.got = run_baseline()
        with open(GOLDEN) as f:
            cls.exp = json.load(f)

    def _check(self, key):
        errors, stats = [], {"n": 0, "max_abs_diff": 0.0}
        compare(key, self.got[key], self.exp[key], errors, stats)
        self.assertFalse(errors, "\n".join(errors[:40]))
        print(f"  [{key}] {stats['n']} numbers, max |diff| = "
              f"{stats['max_abs_diff']:.3e} (backend {self.got['backend']})")

    def test_fields_small_grid(self):
        self._check("fields")

    def test_checksums_production_grids(self):
        self._check("checks")

    def test_reference_model_fields(self):
        self._check("ref")

    def test_fields_are_finite_and_nontrivial(self):
        h = self.got["fields"]["32x16_t1.7"]["h"]
        self.assertTrue(all(math.isfinite(v) for v in h))
        self.assertGreater(max(abs(v) for v in h), 1e-3)


if __name__ == "__main__":
    if "--update" in sys.argv:
        data = run_baseline()
        data.pop("bench", None)
        with open(GOLDEN, "w") as f:
            json.dump(data, f, indent=1)
        print(f"wrote {GOLDEN}")
    else:
        unittest.main(verbosity=2)
