# FFTW usage audit

Scope: every `fftw_*` call in this repository, as a basis for replacing FFTW
(GPL-2.0-or-later) with PocketFFT (BSD-3-Clause). This document records the
state of the code **before** the port; no code was changed for the audit.

Branch audited: `ament_environment_hooks` at `2d2035f`.

## 1. Summary

* FFTW is used in exactly one component: the FFT wave model
  `LinearRandomFFTWaveSimulation` (production) and its reference twin
  `LinearRandomFFTWaveSimulationRef` (tests only), plus one FFTW API unit
  test (`EigenFFTW_TEST.cc`). No plugin, sampler, physics or rendering code
  calls FFTW directly; they use the wave model through `IWaveSimulation`.
* All transforms are **2-D inverse (backward) DFTs in double precision**,
  unnormalised, on row-major `nx x ny` grids. `nx`, `ny` are the
  `cell_count` wave parameters (default `128 x 128`, documented as powers
  of two). There are no forward transforms, no single-precision transforms
  and no 1-D transforms outside the API unit test.
* Plans are created **once, in the constructor**, with `FFTW_ESTIMATE`, and
  executed per simulation step. The production model executes each plan
  lazily, at most once per `SetTime`; the reference model executes on every
  accessor call.
* The production wave-height transform (and the seven derivative /
  displacement transforms and the pressure transforms) **already use
  complex-to-real** (`fftw_plan_dft_c2r_2d`) on the half spectrum
  `nx x (ny/2 + 1)`. The Hermitian symmetry of the spectrum is exploited when
  the amplitudes are built (`ComputeCurrentAmplitudes` fills only
  `iky <= ny/2`). There is no c2c-for-real-output waste in the production
  path.
* The reference model uses **complex-to-complex** (`fftw_plan_dft_2d`,
  `FFTW_BACKWARD`) on the full `nx x ny` spectrum and discards the imaginary
  part of the output. This is the only c2c-where-c2r-suffices site; it runs
  only in unit tests (see §5 for the recommendation).

## 2. Call-site inventory

### 2.1 `gz-waves/src/LinearRandomFFTWaveSimulation.cc` (+ `LinearRandomFFTWaveSimulationImpl.hh`) - production

| Plan | Input array (complex, row-major) | Output array (real, row-major) | API | Shape | Consumer |
|---|---|---|---|---|---|
| `fft_plan0_` | `fft_h_` `nx x (ny/2+1)` | `fft_out0_` `nx x ny` | `fftw_plan_dft_c2r_2d(nx, ny, ...)` | 2-D c2r, double | `ElevationAt` (array and scalar) |
| `fft_plan1_` | `fft_h_ikx_` | `fft_out1_` | c2r 2-D | idem | `ElevationDerivAt` (assigned to `dhdy`) |
| `fft_plan2_` | `fft_h_iky_` | `fft_out2_` | c2r 2-D | idem | `ElevationDerivAt` (assigned to `dhdx`) |
| `fft_plan3_` | `fft_sx_` | `fft_out3_` | c2r 2-D | idem | `DisplacementAt` (`sy = -lambda * out3`) |
| `fft_plan4_` | `fft_sy_` | `fft_out4_` | c2r 2-D | idem | `DisplacementAt` (`sx = -lambda * out4`) |
| `fft_plan5_` | `fft_h_kxkx_` | `fft_out5_` | c2r 2-D | idem | `DisplacementDerivAt` (`dsydy`) |
| `fft_plan6_` | `fft_h_kyky_` | `fft_out6_` | c2r 2-D | idem | `DisplacementDerivAt` (`dsxdx`) |
| `fft_plan7_` | `fft_h_kxky_` | `fft_out7_` | c2r 2-D | idem | `DisplacementDerivAt` (`dsxdy`) |
| `fft_plan_p_[iz]`, `iz < nz` | `fft_in_p_[iz]` | `fft_out_p_[iz]` | c2r 2-D | idem | `PressureAt(iz, ...)` and scalar `PressureAt` |

(The x/y naming swap between input and consumer follows the row/column
convention of the model, where `ikx` indexes rows; it is not an FFT matter.)

Details:

* **Plan lifetime.** `CreateFFTWPlans()` is called from both constructors
  and only there. `SetLambda`, `SetWindVelocity` and `ComputeBaseAmplitudes`
  recompute amplitudes but keep the plans. Flags: `FFTW_ESTIMATE` (no
  measurement, deterministic plan choice).
* **Execution per step.** `SetTime(t)` calls `ComputeCurrentAmplitudes(t)`,
  which rewrites all input arrays and sets every entry of
  `fft_needs_update_` (size `8 + nz`) to true. Each accessor runs its
  transform only if its flag is set, so per step at most 8 + nz transforms
  are executed, and only those whose outputs are requested. The Gazebo
  plugins request elevation, displacements and all derivatives every step
  (`OceanTile::Update` -> `DisplacementAndDerivAt`): 8 transforms of
  `nx x ny` per tile per step, pressure transforms only when `nz > 1`.
* **Threading.** With `libfftw3_omp` found, CMake defines `USE_FFTW3_OMP` and
  the constructor calls `fftw_init_threads()` and
  `fftw_plan_with_nthreads(omp_get_max_threads())`. Without it the
  transforms are single-threaded.
* **Destruction.** `DestroyFFTWPlans()` destroys the eight fixed plans. The
  `nz` pressure plans in `fft_plan_p_` are **never destroyed** (leak, one
  set per model instance). Under `USE_FFTW3_OMP` the destructor also calls
  `fftw_cleanup_threads()`, which is process-global: if another model
  instance is still alive (the waves visual and the ocean tile each own
  one), its plans become invalid. Neither problem has been observed in
  practice because instances are long-lived, but both disappear with a
  backend that has no global state.
* **Input destruction.** FFTW multi-dimensional out-of-place c2r plans
  destroy their input array on execution (`FFTW_PRESERVE_INPUT` is not
  supported for them). The code copes through the lazy flags, so each input
  is transformed exactly once after it is written; the `\todo` in the scalar
  `ElevationAt(ix, iy, eta)` refers to this.
* **Normalisation.** Unnormalised backward transform; the `1/(nx ny)` factor
  is folded into the amplitudes (`cap_psi_norm`, `delta_kx * delta_ky`).
* **Nyquist handling.** The `ikx == nx/2` and `iky == ny/2` derivative terms
  are zeroed before the transform (`fft-deriv.pdf` convention), which keeps
  the c2r input consistent with a real output.

### 2.2 `gz-waves/src/LinearRandomFFTWaveSimulationRef.cc` (+ `...RefImpl.hh`) - reference, tests only

| Plan | Input (complex `nx x ny`) | Output (complex `nx x ny`) | API |
|---|---|---|---|
| `fft_plan0_` .. `fft_plan7_` | `fft_h_`, `fft_h_ikx_`, `fft_h_iky_`, `fft_sx_`, `fft_sy_`, `fft_h_kxkx_`, `fft_h_kyky_`, `fft_h_kxky_` | `fft_out0_` .. `fft_out7_` | `fftw_plan_dft_2d(nx, ny, in, out, FFTW_BACKWARD, FFTW_ESTIMATE)`, 2-D c2c, double |

* Plans created once in the constructor (`CreateFFTWPlans()` after
  `ComputeBaseAmplitudes()`); destroyed in the destructor (all eight).
* Executed **unconditionally on every accessor call** (no lazy flags); the
  accessors take `.real()` of the output. The imaginary part, which is zero
  up to rounding for a Hermitian spectrum, is computed and discarded.
* Used by `LinearRandomFFTWaveSimulation_TEST.cc` (Hermitian-symmetry checks
  that read the full `nx x ny` spectrum arrays, Parseval checks, lambda = 0
  check), and by `WaveSpectrum_TEST.cc` / `WaveSpreadingFunction_TEST.cc`
  for its static spectrum functions only (no FFT). It is compiled into
  `libgz-waves1` but no plugin or sampler instantiates it.
* No threading calls.

### 2.3 `gz-waves/src/EigenFFTW_TEST.cc` - FFTW API unit test

| Test | API | Shape | Notes |
|---|---|---|---|
| `DFT_C2C_1D` | `fftw_plan_dft_1d(n, in, out, FFTW_BACKWARD, FFTW_ESTIMATE)`, three storage variants (`fftw_malloc`, `std::vector`, `Eigen::ArrayXcd`) | 1-D c2c, n = 8, double | expected values from `numpy.fft` with `norm="forward"`; checks real part and zero imaginary part to 1e-15 |
| `DFT_C2R_1D` | `fftw_plan_dft_c2r_1d(n, in, out, FFTW_ESTIMATE)` | 1-D c2r, n = 8 (input n/2+1 = 5), double | same data |
| `EigenStorageOrdering` | none (`#if 0`) | - | Eigen storage demo |

### 2.4 Build system and packaging

| File | FFTW reference |
|---|---|
| `gz-waves/CMakeLists.txt` 318-351 | `find_path(FFTW3-3_INCLUDE_DIR)`, `find_library(FFTW3-3_LIBRARY)`, optional `find_library(FFTW3_OMP_LIBRARY)`; `add_definitions(-DUSE_FFTW3)`, `-DUSE_FFTW3_OMP`; `FFT_INCLUDE_DIRS`, `FFT_LIBRARIES` |
| `gz-waves/src/CMakeLists.txt` 65, 93 | `${FFT_LIBRARIES}` linked `PUBLIC` into `gz-waves1` and into the unit tests |
| `gz-waves/src/CMakeLists.txt` 40 | `EigenFFTW_TEST.cc` in `gtest_sources` |
| `README.md` 21, 32-36, 54-60, 68-72 | licence statement (GPL), `libfftw3-dev`, `brew install ... fftw`, threaded-FFTW note |
| `.github/workflows/ubuntu-jammy-ci.yml` 37 | `apt-get install ... libfftw3-dev` |
| `.github/scripts/brew_install_deps.sh` 27 | `fftw` |
| `LICENSE_THIRDPARTY` | no FFTW entry (the file lists Boost, Eigen and CGAL texts only) |
| Docker | **no Dockerfile or compose file exists in this repository**; nothing to change |

`fftw3.h` is included by the two `Impl.hh` headers, which live in `src/` and
are not installed, so FFTW does not appear in the public API. `<omp.h>` is
included in `LinearRandomFFTWaveSimulation.cc` only under `USE_FFTW3_OMP`.

## 3. Transform shapes actually exercised

| Context | `nx x ny` | Transforms per step | Precision |
|---|---|---|---|
| Gazebo waves visual / ocean tile, default `cell_count` | 128 x 128 (`tile_size` 256 m) | 8 c2r (+ `nz` pressure c2r when `nz > 1`) | double |
| Unit tests (`LinearRandomFFTWaveSimFixture`) | small powers of two (4 x 4 .. 32 x 32) | up to 8 | double |
| Performance tests (`test/performance`) | 128 x 128, 256 x 256 | 8 | double |

No caller uses `float`; the float tolerance of the validation plan is
therefore not exercised by this repository.

## 4. Hermitian symmetry and c2c vs c2r

The elevation spectrum `zhat_` is constructed so that
`zhat(-kx, -ky) = conj(zhat(kx, ky))` (the loops in `ComputeCurrentAmplitudes`
set the conjugate pairs explicitly and zero the DC term). A backward DFT of a
Hermitian spectrum is real, so a complex-to-real transform of the half
spectrum `iky = 0 .. ny/2` reproduces the full complex-to-complex result's
real part exactly (up to rounding) at half the arithmetic and half the input
storage.

* **Production model:** already c2r for all 8 + nz transforms. No change
  available.
* **Reference model:** c2c on the full spectrum, imaginary part discarded.
  A c2r transform would halve its FFT work. Its full `nx x ny` spectrum
  arrays are, however, what the eight `Hermitian*Reference` unit tests
  inspect, and the reference's purpose is to be the straightforward,
  unoptimised implementation that the optimised class is checked against.

## 5. Recommendation for the port

1. Introduce an internal `fft::` interface (backward c2r and c2c transforms
   on row-major arrays of a given shape, plan held by the object, input
   preserved) and move both models and the API unit test onto it with the
   FFTW backend, verifying no behaviour change.
2. Swap the backend to PocketFFT (`pocketfft_hdronly.h`, BSD-3-Clause,
   header-only): same sign convention (`forward=false` is `e^{+i...}`), same
   unnormalised scaling with `fct = 1`, same half-spectrum layout along the
   last axis, input `const` (not destroyed). Enable its plan cache so the
   twiddle tables are computed once per length, as FFTW plans were.
3. Keep the reference model complex-to-complex: it is test-only and its
   value is independence from the optimised path. Report this explicitly in
   the PR rather than switching it silently (the alternative, c2r on a
   strided half-plane view of the same arrays, is straightforward if
   preferred).
4. Fix the two lifecycle issues in passing, since the interface owns its
   plans: pressure plans are destroyed with the object, and no
   process-global thread clean-up is needed.
5. Remove FFTW from CMake, README, CI and the brew script; add the PocketFFT
   licence to `LICENSE_THIRDPARTY`.
