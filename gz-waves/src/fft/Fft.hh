// Copyright (C) 2026  Ethan Chiang and contributors
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

/// \file fft/Fft.hh
/// \brief Internal interface for the discrete Fourier transforms used by
/// the FFT wave models. Hides the FFT library behind two plan objects so
/// the backend can be swapped without touching the models
/// (docs/fftw_audit.md).
///
/// Conventions (identical to FFTW's backward transform and to
/// numpy.fft.ifft with norm="forward"):
///   out[j] = sum_k in[k] * exp(+2 pi i j.k / n)   (no 1/n factor)
/// Arrays are contiguous and row-major; `shape` lists the dimensions
/// outermost first. For the complex-to-real transform the input has the
/// last dimension halved, shape[..., n_last / 2 + 1] (the non-negative
/// frequencies), and the output has the full `shape`.

#ifndef GZ_WAVES_SRC_FFT_FFT_HH_
#define GZ_WAVES_SRC_FFT_FFT_HH_

#include <complex>
#include <memory>
#include <vector>

#include "gz/waves/Types.hh"

namespace gz
{
namespace waves
{
namespace fft
{
/// \brief Shape of a row-major array, outermost dimension first.
using Shape = std::vector<Index>;

/// \brief Name of the FFT backend compiled into the library.
const char* Backend();

/// \brief Backward complex-to-real transform bound to an input/output pair.
///
/// The plan is created once, in the constructor, and Execute() may be
/// called any number of times. Both arrays must outlive the plan.
class BackwardC2R
{
 public:
  /// \param[in] _shape Shape of the real output array.
  /// \param[in] _in    Complex input, shape[..., n_last / 2 + 1].
  /// \param[out] _out  Real output, `_shape`.
  BackwardC2R(const Shape& _shape, const std::complex<double>* _in,
      double* _out);
  ~BackwardC2R();
  BackwardC2R(BackwardC2R&& _other) noexcept;
  BackwardC2R& operator=(BackwardC2R&& _other) noexcept;
  BackwardC2R(const BackwardC2R&) = delete;
  BackwardC2R& operator=(const BackwardC2R&) = delete;

  /// \brief Run the transform: read the input array, write the output.
  /// \note Whether the input survives Execute() depends on the backend;
  /// callers must rewrite it before the next Execute().
  void Execute();

 private:
  class Impl;
  std::unique_ptr<Impl> impl_;
};

/// \brief Backward complex-to-complex transform bound to an input/output
/// pair. Same lifetime rules as BackwardC2R.
class BackwardC2C
{
 public:
  /// \param[in] _shape Shape of both arrays.
  BackwardC2C(const Shape& _shape, const std::complex<double>* _in,
      std::complex<double>* _out);
  ~BackwardC2C();
  BackwardC2C(BackwardC2C&& _other) noexcept;
  BackwardC2C& operator=(BackwardC2C&& _other) noexcept;
  BackwardC2C(const BackwardC2C&) = delete;
  BackwardC2C& operator=(const BackwardC2C&) = delete;

  /// \brief Run the transform.
  void Execute();

 private:
  class Impl;
  std::unique_ptr<Impl> impl_;
};

}  // namespace fft
}  // namespace waves
}  // namespace gz

#endif  // GZ_WAVES_SRC_FFT_FFT_HH_
