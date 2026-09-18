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
/// \brief FFTW3 backend for fft::BackwardC2R / fft::BackwardC2C.

#include "fft/Fft.hh"

#include <fftw3.h>

#ifdef USE_FFTW3_OMP
#include <omp.h>
#endif

#include <complex>
#include <vector>

namespace gz
{
namespace waves
{
namespace fft
{
namespace
{
/// \brief Configure FFTW threading once per process (idempotent).
void InitThreadsOnce()
{
#ifdef USE_FFTW3_OMP
  static const bool done = []()
  {
    fftw_init_threads();
    fftw_plan_with_nthreads(omp_get_max_threads());
    return true;
  }();
  (void)done;
#endif
}

std::vector<int> ToInt(const Shape& _shape)
{
  return std::vector<int>(_shape.begin(), _shape.end());
}
}  // namespace

//////////////////////////////////////////////////
const char* Backend()
{
  return "FFTW3";
}

//////////////////////////////////////////////////
class BackwardC2R::Impl
{
 public:
  Impl(const Shape& _shape, const std::complex<double>* _in, double* _out)
  {
    InitThreadsOnce();
    const std::vector<int> n = ToInt(_shape);
    // FFTW may overwrite the input of an out-of-place c2r plan; the
    // pointer is only const at the interface.
    plan = fftw_plan_dft_c2r(static_cast<int>(n.size()), n.data(),
        reinterpret_cast<fftw_complex*>(
            const_cast<std::complex<double>*>(_in)),
        _out, FFTW_ESTIMATE);
  }
  ~Impl()
  {
    if (plan)
      fftw_destroy_plan(plan);
  }
  fftw_plan plan{nullptr};
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
  fftw_execute(impl_->plan);
}

//////////////////////////////////////////////////
class BackwardC2C::Impl
{
 public:
  Impl(const Shape& _shape, const std::complex<double>* _in,
      std::complex<double>* _out)
  {
    InitThreadsOnce();
    const std::vector<int> n = ToInt(_shape);
    plan = fftw_plan_dft(static_cast<int>(n.size()), n.data(),
        reinterpret_cast<fftw_complex*>(
            const_cast<std::complex<double>*>(_in)),
        reinterpret_cast<fftw_complex*>(_out),
        FFTW_BACKWARD, FFTW_ESTIMATE);
  }
  ~Impl()
  {
    if (plan)
      fftw_destroy_plan(plan);
  }
  fftw_plan plan{nullptr};
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
  fftw_execute(impl_->plan);
}

}  // namespace fft
}  // namespace waves
}  // namespace gz
