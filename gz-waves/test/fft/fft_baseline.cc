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

/// \file fft_baseline.cc
/// \brief Golden-value tool for the FFT backend replacement
/// (docs/fftw_audit.md). Prints JSON on stdout:
///
///   - "fields":  full elevation / displacement / derivative / pressure
///                fields of the FFT wave model for a small grid at fixed
///                spectrum parameters and times (deterministic: the model
///                seeds its random amplitudes with a fixed seed).
///   - "checks":  sum, sum of squares, max |.| and 64 sampled values of
///                the same fields on the production grid sizes.
///   - "ref":     full fields of the reference (c2c) model.
///   - "bench":   with --bench N, per-step timing in microseconds.
///
/// Usage: fft_baseline [--bench N]

#include <Eigen/Dense>

#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <functional>
#include <memory>
#include <string>
#include <vector>

#include <gz/waves/LinearRandomFFTWaveSimulation.hh>
#include <gz/waves/Types.hh>

#include "LinearRandomFFTWaveSimulationRefImpl.hh"
#include "fft/Fft.hh"

using gz::waves::Index;
using gz::waves::LinearRandomFFTWaveSimulation;
using gz::waves::LinearRandomFFTWaveSimulationRef;

namespace
{
//////////////////////////////////////////////////
// Minimal JSON writer (numbers with 17 significant digits).
class Json
{
 public:
  void Begin(const char* key = nullptr) { Key(key); Put("{"); first_ = true; }
  void End() { Put("}"); first_ = false; }
  void Num(const char* key, double v)
  {
    Key(key);
    char buf[64];
    if (std::isfinite(v))
      std::snprintf(buf, sizeof(buf), "%.17g", v);
    else
      std::snprintf(buf, sizeof(buf), "null");
    Put(buf);
    first_ = false;
  }
  void Str(const char* key, const std::string& v)
  {
    Key(key); Put("\"" + v + "\""); first_ = false;
  }
  void Arr(const char* key, const Eigen::ArrayXXd& a)
  {
    Key(key);
    std::string s = "[";
    char buf[64];
    for (Index i = 0; i < a.size(); ++i)
    {
      std::snprintf(buf, sizeof(buf), "%.17g", a(i));
      if (i) s += ",";
      s += buf;
    }
    s += "]";
    Put(s);
    first_ = false;
  }
  void Arr(const char* key, const std::vector<double>& a)
  {
    Eigen::ArrayXXd e(a.size(), 1);
    for (size_t i = 0; i < a.size(); ++i) e(i, 0) = a[i];
    Arr(key, e);
  }

 private:
  void Key(const char* key)
  {
    if (!first_) std::fputs(",\n", stdout);
    if (key) std::printf("\"%s\": ", key);
  }
  void Put(const std::string& s) { std::fputs(s.c_str(), stdout); }
  bool first_{true};
};

//////////////////////////////////////////////////
struct Fields
{
  Eigen::ArrayXXd h, sx, sy, dhdx, dhdy, dsxdx, dsydy, dsxdy;
  std::vector<Eigen::ArrayXXd> pressure;
};

Fields Compute(const LinearRandomFFTWaveSimulation& sim, Index nz)
{
  const Index n2 = sim.SizeX() * sim.SizeY();
  Fields f;
  for (Eigen::ArrayXXd* a : {&f.h, &f.sx, &f.sy, &f.dhdx, &f.dhdy,
      &f.dsxdx, &f.dsydy, &f.dsxdy})
    *a = Eigen::ArrayXXd::Zero(n2, 1);
  sim.DisplacementAndDerivAt(f.h, f.sx, f.sy, f.dhdx, f.dhdy,
      f.dsxdx, f.dsydy, f.dsxdy);
  for (Index iz = 0; iz < nz; ++iz)
  {
    Eigen::ArrayXXd p = Eigen::ArrayXXd::Zero(n2, 1);
    sim.PressureAt(iz, p);
    f.pressure.push_back(p);
  }
  return f;
}

void WriteFields(Json& j, const Fields& f)
{
  j.Arr("h", f.h);
  j.Arr("sx", f.sx);
  j.Arr("sy", f.sy);
  j.Arr("dhdx", f.dhdx);
  j.Arr("dhdy", f.dhdy);
  j.Arr("dsxdx", f.dsxdx);
  j.Arr("dsydy", f.dsydy);
  j.Arr("dsxdy", f.dsxdy);
  for (size_t iz = 0; iz < f.pressure.size(); ++iz)
    j.Arr(("pressure_z" + std::to_string(iz)).c_str(), f.pressure[iz]);
}

void WriteChecks(Json& j, const char* key, const Eigen::ArrayXXd& a)
{
  j.Begin(key);
  j.Num("sum", a.sum());
  j.Num("sum_sq", (a * a).sum());
  j.Num("max_abs", a.abs().maxCoeff());
  // 64 samples spread over the field (fixed stride, so the same grid
  // points are compared across backends).
  std::vector<double> s;
  const Index stride = std::max<Index>(1, a.size() / 64);
  for (Index i = 0; i < a.size() && s.size() < 64; i += stride)
    s.push_back(a(i));
  j.Arr("samples", s);
  j.End();
}

double Now()
{
  using clock = std::chrono::steady_clock;
  return std::chrono::duration<double, std::micro>(
      clock::now().time_since_epoch()).count();
}
}  // namespace

//////////////////////////////////////////////////
int main(int argc, char** argv)
{
  int bench = 0;
  for (int i = 1; i < argc; ++i)
  {
    if (std::strcmp(argv[i], "--bench") == 0 && i + 1 < argc)
      bench = std::atoi(argv[++i]);
  }

  // Common spectrum parameters (tile 256 m, wind 8 m/s at 30 deg,
  // lambda 0.6, three pressure levels down to 10 m).
  const double lx = 256.0, ly = 256.0, lz = 10.0;
  const double ux = 8.0 * std::cos(M_PI / 6.0);
  const double uy = 8.0 * std::sin(M_PI / 6.0);
  const double lambda = 0.6;
  const std::vector<double> times{0.0, 1.7, 12.3};

  Json j;
  j.Begin();
  j.Str("backend", gz::waves::fft::Backend());

  // Full fields on a small grid.
  j.Begin("fields");
  {
    const Index nx = 32, ny = 16, nz = 3;
    LinearRandomFFTWaveSimulation sim(lx, ly, lz, nx, ny, nz);
    sim.SetWindVelocity(ux, uy);
    sim.SetLambda(lambda);
    for (double t : times)
    {
      sim.SetTime(t);
      Fields f = Compute(sim, nz);
      char key[64];
      std::snprintf(key, sizeof(key), "%ldx%ld_t%.1f",
          static_cast<int64_t>(nx), static_cast<int64_t>(ny), t);
      j.Begin(key);
      WriteFields(j, f);
      j.End();
    }
  }
  j.End();

  // Checksums on production-size grids.
  j.Begin("checks");
  for (Index n : {128, 256})
  {
    LinearRandomFFTWaveSimulation sim(lx, ly, n, n);
    sim.SetWindVelocity(ux, uy);
    sim.SetLambda(lambda);
    sim.SetTime(1.7);
    Fields f = Compute(sim, 0);
    char key[64];
    std::snprintf(key, sizeof(key), "%ldx%ld_t1.7",
        static_cast<int64_t>(n), static_cast<int64_t>(n));
    j.Begin(key);
    WriteChecks(j, "h", f.h);
    WriteChecks(j, "sx", f.sx);
    WriteChecks(j, "sy", f.sy);
    WriteChecks(j, "dhdx", f.dhdx);
    WriteChecks(j, "dhdy", f.dhdy);
    WriteChecks(j, "dsxdx", f.dsxdx);
    WriteChecks(j, "dsydy", f.dsydy);
    WriteChecks(j, "dsxdy", f.dsxdy);
    j.End();
  }
  j.End();

  // Reference (complex-to-complex) model, full fields.
  j.Begin("ref");
  {
    const Index nx = 16, ny = 8;
    // The public wrapper does not implement the full IWaveSimulation
    // interface; use the implementation class as the unit tests do.
    LinearRandomFFTWaveSimulationRef::Impl ref(200.0, 100.0, nx, ny);
    ref.SetWindVelocity(ux, uy);
    ref.SetLambda(lambda);
    for (double t : {0.0, 1.7})
    {
      ref.SetTime(t);
      const Index n2 = nx * ny;
      Fields f;
      for (Eigen::ArrayXXd* a : {&f.h, &f.sx, &f.sy, &f.dhdx, &f.dhdy,
          &f.dsxdx, &f.dsydy, &f.dsxdy})
        *a = Eigen::ArrayXXd::Zero(n2, 1);
      ref.ElevationAt(f.h);
      ref.ElevationDerivAt(f.dhdx, f.dhdy);
      ref.DisplacementAt(f.sx, f.sy);
      ref.DisplacementDerivAt(f.dsxdx, f.dsydy, f.dsxdy);
      char key[64];
      std::snprintf(key, sizeof(key), "%ldx%ld_t%.1f",
          static_cast<int64_t>(nx), static_cast<int64_t>(ny), t);
      j.Begin(key);
      WriteFields(j, f);
      j.End();
    }
  }
  j.End();

  // Timing: per simulation step, amplitude update (SetTime) and the eight
  // backward transforms plus copies (DisplacementAndDerivAt), and the
  // eight transforms alone (a second call after SetTime re-runs nothing,
  // so the transform cost is the difference between a call with dirty
  // flags and a call without).
  if (bench > 0)
  {
    j.Begin("bench");
    for (Index n : {128, 256, 512})
    {
      LinearRandomFFTWaveSimulation sim(lx, ly, n, n);
      sim.SetWindVelocity(ux, uy);
      sim.SetLambda(lambda);
      const Index n2 = n * n;
      Fields f;
      for (Eigen::ArrayXXd* a : {&f.h, &f.sx, &f.sy, &f.dhdx, &f.dhdy,
          &f.dsxdx, &f.dsydy, &f.dsxdy})
        *a = Eigen::ArrayXXd::Zero(n2, 1);
      auto all = [&]()
      {
        sim.DisplacementAndDerivAt(f.h, f.sx, f.sy, f.dhdx, f.dhdy,
            f.dsxdx, f.dsydy, f.dsxdy);
      };
      // warm up
      sim.SetTime(0.5); all();

      double tAmp = 0.0, tFft = 0.0, tCopy = 0.0;
      for (int i = 0; i < bench; ++i)
      {
        const double t = 0.01 * i;
        double t0 = Now();
        sim.SetTime(t);
        double t1 = Now();
        all();              // 8 transforms + copies
        double t2 = Now();
        all();              // copies only (flags clean)
        double t3 = Now();
        tAmp += t1 - t0;
        tFft += (t2 - t1) - (t3 - t2);
        tCopy += t3 - t2;
      }
      char key[64];
      std::snprintf(key, sizeof(key), "%ldx%ld",
          static_cast<int64_t>(n), static_cast<int64_t>(n));
      j.Begin(key);
      j.Num("steps", bench);
      j.Num("amplitudes_us_per_step", tAmp / bench);
      j.Num("fft8_us_per_step", tFft / bench);
      j.Num("fft_us_per_transform", tFft / bench / 8.0);
      j.Num("copy_us_per_step", tCopy / bench);
      j.End();
    }
    j.End();
  }

  j.End();
  std::printf("\n");
  return 0;
}
