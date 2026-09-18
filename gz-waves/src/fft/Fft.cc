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

/// \file fft/Fft.cc
/// \brief PocketFFT backend for fft::BackwardC2R / fft::BackwardC2C.
///
/// PocketFFT (thirdparty/pocketfft, BSD-3-Clause) is header-only. Its
/// c2r / c2c functions take the array shape, byte strides and the axes to
/// transform; `forward = false` is the e^{+i...} (backward) transform and
/// `fct = 1` leaves it unnormalised, matching the previous FFTW3 plans.
/// The plan cache (POCKETFFT_CACHE_SIZE) keeps the twiddle tables of the
/// most recent transform lengths alive across calls, so repeated
/// transforms of the same size do no planning work, as FFTW plans did not.
/// PocketFFT does not modify its input.

#include "fft/Fft.hh"

#include <complex>
#include <cstddef>
#include <vector>

#define POCKETFFT_CACHE_SIZE 16
#include "pocketfft_hdronly.h"

namespace gz
{
namespace waves
{
namespace fft
{
namespace
{
/// \brief Threads per transform. PocketFFT can split one multi-dimensional
/// transform over a std::thread pool, but for the grids used here
/// (128^2 .. 512^2) the pool overhead exceeds the gain: on 4 cores a
/// 128x128 c2r took 139 us with 4 threads against 66 us with one. The
/// wave model instead runs its eight independent per-step transforms
/// concurrently (LinearRandomFFTWaveSimulation::Impl::ExecutePending), so
/// each transform stays single-threaded.
constexpr size_t kThreads = 1;

/// \brief Row-major byte strides of a contiguous array of `_shape` whose
/// elements are `_elem` bytes.
pocketfft::stride_t ContiguousStrides(const pocketfft::shape_t& _shape,
    size_t _elem)
{
  pocketfft::stride_t s(_shape.size());
  ptrdiff_t stride = static_cast<ptrdiff_t>(_elem);
  for (size_t i = _shape.size(); i-- > 0;)
  {
    s[i] = stride;
    stride *= static_cast<ptrdiff_t>(_shape[i]);
  }
  return s;
}

pocketfft::shape_t ToShape(const Shape& _shape)
{
  return pocketfft::shape_t(_shape.begin(), _shape.end());
}

pocketfft::shape_t AllAxes(size_t _rank)
{
  pocketfft::shape_t axes(_rank);
  for (size_t i = 0; i < _rank; ++i)
    axes[i] = i;
  return axes;
}
}  // namespace

//////////////////////////////////////////////////
const char* Backend()
{
  return "PocketFFT";
}

//////////////////////////////////////////////////
class BackwardC2R::Impl
{
 public:
  Impl(const Shape& _shape, const std::complex<double>* _in, double* _out) :
    shapeOut(ToShape(_shape)),
    axes(AllAxes(_shape.size())),
    in(_in),
    out(_out)
  {
    pocketfft::shape_t shapeIn = shapeOut;
    shapeIn.back() = shapeOut.back() / 2 + 1;
    strideIn = ContiguousStrides(shapeIn, sizeof(std::complex<double>));
    strideOut = ContiguousStrides(shapeOut, sizeof(double));
  }

  void Execute()
  {
    pocketfft::c2r<double>(shapeOut, strideIn, strideOut, axes,
        /*forward=*/false, in, out, /*fct=*/1.0, kThreads);
  }

  pocketfft::shape_t shapeOut;
  pocketfft::shape_t axes;
  pocketfft::stride_t strideIn;
  pocketfft::stride_t strideOut;
  const std::complex<double>* in;
  double* out;
};

BackwardC2R::BackwardC2R(const Shape& _shape,
    const std::complex<double>* _in, double* _out) :
  impl_(new Impl(_shape, _in, _out))
{
}
BackwardC2R::~BackwardC2R() = default;
BackwardC2R::BackwardC2R(BackwardC2R&&) noexcept = default;
BackwardC2R& BackwardC2R::operator=(BackwardC2R&&) noexcept = default;

void BackwardC2R::Execute()
{
  impl_->Execute();
}

//////////////////////////////////////////////////
class BackwardC2C::Impl
{
 public:
  Impl(const Shape& _shape, const std::complex<double>* _in,
      std::complex<double>* _out) :
    shape(ToShape(_shape)),
    axes(AllAxes(_shape.size())),
    stride(ContiguousStrides(shape, sizeof(std::complex<double>))),
    in(_in),
    out(_out)
  {
  }

  void Execute()
  {
    pocketfft::c2c<double>(shape, stride, stride, axes,
        /*forward=*/false, in, out, /*fct=*/1.0, kThreads);
  }

  pocketfft::shape_t shape;
  pocketfft::shape_t axes;
  pocketfft::stride_t stride;
  const std::complex<double>* in;
  std::complex<double>* out;
};

BackwardC2C::BackwardC2C(const Shape& _shape,
    const std::complex<double>* _in, std::complex<double>* _out) :
  impl_(new Impl(_shape, _in, _out))
{
}
BackwardC2C::~BackwardC2C() = default;
BackwardC2C::BackwardC2C(BackwardC2C&&) noexcept = default;
BackwardC2C& BackwardC2C::operator=(BackwardC2C&&) noexcept = default;

void BackwardC2C::Execute()
{
  impl_->Execute();
}

}  // namespace fft
}  // namespace waves
}  // namespace gz
