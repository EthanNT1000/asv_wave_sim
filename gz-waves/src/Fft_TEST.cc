// Copyright (C) 2022  Rhys Mainwaring
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

/// \file Fft_TEST.cc
/// \brief Checks of the fft:: interface against values from numpy.fft
/// (independent of the backend compiled in).

#include <gtest/gtest.h>

#include <Eigen/Dense>

#include <complex>
#include <string>
#include <vector>

#include "fft/Fft.hh"

using gz::waves::Index;
namespace fft = gz::waves::fft;

namespace
{
// Python code to generate test data
//
// x = [0.0, 0.0, 0.0, 1.0, 1.0, 0.0, 0.0, 0.0]
// xhat = fft.fft(x, 8, norm="forward")
// xx = fft.ifft(xhat, 8, norm="forward")
//
// with np.printoptions(precision=16, suppress=True):
//     print(f"x:\n{x}")
//     print(f"xhat:\n{xhat.real}\n{xhat.imag}")
//     print(f"xx:\n{xx.real}\n{xx.imag}")
//
constexpr Index kN = 8;

Eigen::ArrayXcd ExpectedX()
{
  Eigen::ArrayXcd x = Eigen::ArrayXcd::Zero(kN);
  x.real() << 0.0, 0.0, 0.0, 1.0, 1.0, 0.0, 0.0, 0.0;
  return x;
}

Eigen::ArrayXcd ExpectedXHat()
{
  Eigen::ArrayXcd xhat = Eigen::ArrayXcd::Zero(kN);
  xhat.real() <<  0.25,   -0.2133883476483184,  0.125,  -0.0366116523516816,
                  0.,     -0.0366116523516816,  0.125,  -0.2133883476483184;
  xhat.imag() <<  0.,     -0.0883883476483184,  0.125,  -0.0883883476483184,
                  0.,      0.0883883476483184, -0.125,   0.0883883476483184;
  return xhat;
}
}  // namespace

//////////////////////////////////////////////////
TEST(Fft, BackendName)
{
  EXPECT_NE(fft::Backend(), nullptr);
  EXPECT_GT(std::string(fft::Backend()).size(), 0u);
}

//////////////////////////////////////////////////
TEST(Fft, BackwardC2C_1D)
{
  const Eigen::ArrayXcd x = ExpectedX();
  const Eigen::ArrayXcd xhat = ExpectedXHat();

  {  // std::vector storage
    std::vector<std::complex<double>> in(kN);
    std::vector<std::complex<double>> out(kN);
    fft::BackwardC2C plan({kN}, in.data(), out.data());
    for (Index i = 0; i < kN; ++i)
      in[i] = xhat(i);
    plan.Execute();
    for (Index i = 0; i < kN; ++i)
    {
      EXPECT_NEAR(out[i].real(), x(i).real(), 1.0E-15);
      EXPECT_NEAR(out[i].imag(), 0.0, 1.0E-15);
    }
  }

  {  // Eigen storage, executed twice (plan reuse)
    Eigen::ArrayXcd in = Eigen::ArrayXcd::Zero(kN);
    Eigen::ArrayXcd out = Eigen::ArrayXcd::Zero(kN);
    fft::BackwardC2C plan({kN}, in.data(), out.data());
    for (int pass = 0; pass < 2; ++pass)
    {
      in = xhat;
      plan.Execute();
      for (Index i = 0; i < kN; ++i)
      {
        EXPECT_NEAR(out(i).real(), x(i).real(), 1.0E-15);
        EXPECT_NEAR(out(i).imag(), 0.0, 1.0E-15);
      }
    }
  }
}

//////////////////////////////////////////////////
TEST(Fft, BackwardC2R_1D)
{
  const Eigen::ArrayXcd x = ExpectedX();
  const Eigen::ArrayXcd xhat = ExpectedXHat();

  std::vector<std::complex<double>> in(kN/2 + 1);
  std::vector<double> out(kN, 0.0);
  fft::BackwardC2R plan({kN}, in.data(), out.data());
  for (int pass = 0; pass < 2; ++pass)
  {
    for (Index i = 0; i < kN/2 + 1; ++i)
      in[i] = xhat(i);
    plan.Execute();
    for (Index i = 0; i < kN; ++i)
      EXPECT_NEAR(out[i], x(i).real(), 1.0E-15);
  }
}

//////////////////////////////////////////////////
/// 2-D transforms: the c2r result of the half spectrum must equal the real
/// part of the c2c result of the full Hermitian spectrum, and the 1-D
/// results must appear along each row when the other dimension is 1.
TEST(Fft, BackwardC2R_2D_MatchesC2C)
{
  const Index nx = 4, ny = 8;
  using RowMajorC = Eigen::Array<std::complex<double>, Eigen::Dynamic,
      Eigen::Dynamic, Eigen::RowMajor>;
  using RowMajorD = Eigen::Array<double, Eigen::Dynamic,
      Eigen::Dynamic, Eigen::RowMajor>;

  // Hermitian spectrum: z(-kx, -ky) = conj(z(kx, ky)), zero DC.
  RowMajorC full = RowMajorC::Zero(nx, ny);
  for (Index ikx = 0; ikx < nx; ++ikx)
  {
    for (Index iky = 0; iky < ny; ++iky)
    {
      const Index ckx = (nx - ikx) % nx;
      const Index cky = (ny - iky) % ny;
      const bool self = (ckx == ikx && cky == iky);
      const double re = 0.1 * (ikx + 1) + 0.01 * (iky + 1);
      const double im = self ? 0.0 : 0.05 * (ikx - iky);
      full(ikx, iky) = std::complex<double>(re, im);
      full(ckx, cky) = std::complex<double>(re, -im);
    }
  }
  full(0, 0) = 0.0;

  RowMajorC half = full.leftCols(ny/2 + 1);
  RowMajorD outR = RowMajorD::Zero(nx, ny);
  RowMajorC outC = RowMajorC::Zero(nx, ny);

  fft::BackwardC2R planR({nx, ny}, half.data(), outR.data());
  fft::BackwardC2C planC({nx, ny}, full.data(), outC.data());
  planR.Execute();
  planC.Execute();

  for (Index ikx = 0; ikx < nx; ++ikx)
  {
    for (Index iky = 0; iky < ny; ++iky)
    {
      EXPECT_NEAR(outC(ikx, iky).imag(), 0.0, 1.0E-13);
      EXPECT_NEAR(outR(ikx, iky), outC(ikx, iky).real(), 1.0E-13);
    }
  }

  // Direct evaluation of the DFT definition at one point.
  const Index jx = 1, jy = 3;
  std::complex<double> sum(0.0, 0.0);
  for (Index ikx = 0; ikx < nx; ++ikx)
    for (Index iky = 0; iky < ny; ++iky)
      sum += full(ikx, iky) * std::polar(1.0, 2.0 * M_PI *
          (static_cast<double>(jx * ikx) / nx +
           static_cast<double>(jy * iky) / ny));
  EXPECT_NEAR(outR(jx, jy), sum.real(), 1.0E-13);
}

//////////////////////////////////////////////////
int main(int argc, char **argv)
{
  testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}
