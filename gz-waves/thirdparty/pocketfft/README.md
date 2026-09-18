# PocketFFT (vendored)

`pocketfft_hdronly.h` is the header-only C++ PocketFFT by Martin Reinecke,
used by `gz-waves` for the inverse Fourier transforms of the FFT wave model
(see `src/fft/Fft.cc` and `docs/fftw_audit.md`).

* Upstream: https://github.com/mreineck/pocketfft (branch `cpp`)
* Commit: c90e55b3d529f8efa40ed01a20de22405f45fc65 (2026-06-30)
* Licence: BSD-3-Clause, see `LICENSE.md` in this directory (also listed in
  the repository's `LICENSE_THIRDPARTY`).
* Local modifications: none. Configuration is done from `src/fft/Fft.cc`
  with the `POCKETFFT_CACHE_SIZE` macro before the include.

To update, replace `pocketfft_hdronly.h` and `LICENSE.md` with the upstream
files and record the new commit here.
