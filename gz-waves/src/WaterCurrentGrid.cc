// gz-waves/src/WaterCurrentGrid.cc
#include "gz/waves/WaterCurrentGrid.hh"

#include <cstring>
#include <fstream>
#include <gz/common/Console.hh>

namespace gz::waves {

static constexpr uint32_t MAGIC   = 0x57435247u;
static constexpr uint32_t VERSION = 1u;

bool WaterCurrentGrid::LoadFromFile(const std::string& path)
{
    std::ifstream f(path, std::ios::binary);
    if (!f.is_open()) {
        gzerr << "[WaterCurrentGrid] Cannot open: " << path << "\n";
        return false;
    }

    // ── header ──────────────────────────────────────────────────────────────
    uint32_t magic, version;
    int32_t  nx, ny;
    double   x_min, y_min, csx, csy;

    f.read(reinterpret_cast<char*>(&magic),   4);
    f.read(reinterpret_cast<char*>(&version), 4);
    f.read(reinterpret_cast<char*>(&nx),      4);
    f.read(reinterpret_cast<char*>(&ny),      4);
    f.read(reinterpret_cast<char*>(&x_min),   8);
    f.read(reinterpret_cast<char*>(&y_min),   8);
    f.read(reinterpret_cast<char*>(&csx),     8);
    f.read(reinterpret_cast<char*>(&csy),     8);

    if (magic != MAGIC) {
        gzerr << "[WaterCurrentGrid] Invalid magic number in: " << path << "\n";
        return false;
    }
    if (version != VERSION) {
        gzerr << "[WaterCurrentGrid] Unsupported version " << version << "\n";
        return false;
    }
    if (nx <= 0 || ny <= 0 || csx <= 0.0 || csy <= 0.0) {
        gzerr << "[WaterCurrentGrid] Corrupt header in: " << path << "\n";
        return false;
    }

    // ── data ────────────────────────────────────────────────────────────────
    const size_t n_floats = static_cast<size_t>(nx) * ny * 2;
    data_.resize(n_floats);
    f.read(reinterpret_cast<char*>(data_.data()), n_floats * sizeof(float));

    if (!f) {
        gzerr << "[WaterCurrentGrid] Truncated data in: " << path << "\n";
        return false;
    }

    nx_          = nx;
    ny_          = ny;
    x_min_       = x_min;
    y_min_       = y_min;
    cell_size_x_ = csx;
    cell_size_y_ = csy;
    loaded_      = true;

    gzmsg << "[WaterCurrentGrid] Loaded " << nx_ << "x" << ny_
          << " grid from " << path
          << "  X[" << XMin() << ", " << XMax() << "]"
          << "  Y[" << YMin() << ", " << YMax() << "]\n";
    return true;
}

gz::cgal::Vector3 WaterCurrentGrid::SampleAt(double x, double y) const
{
    if (!loaded_) return {0.0, 0.0, 0.0};

    // Map world coords to fractional grid indices
    double fx = (x - x_min_) / cell_size_x_;
    double fy = (y - y_min_) / cell_size_y_;

    // Out-of-bounds → zero current (no extrapolation)
    if (fx < 0.0 || fy < 0.0 || fx > nx_ - 1 || fy > ny_ - 1)
        return {0.0, 0.0, 0.0};

    // Bilinear interpolation
    int ix0 = static_cast<int>(fx);
    int iy0 = static_cast<int>(fy);
    int ix1 = std::min(ix0 + 1, nx_ - 1);
    int iy1 = std::min(iy0 + 1, ny_ - 1);

    double tx = fx - ix0;   // [0, 1)
    double ty = fy - iy0;

    // Helper: fetch (vx, vy) at grid cell (ix, iy)
    auto cell = [&](int ix, int iy) -> std::pair<double, double> {
        const size_t idx = (static_cast<size_t>(iy) * nx_ + ix) * 2;
        return { data_[idx], data_[idx + 1] };
    };

    auto [vx00, vy00] = cell(ix0, iy0);
    auto [vx10, vy10] = cell(ix1, iy0);
    auto [vx01, vy01] = cell(ix0, iy1);
    auto [vx11, vy11] = cell(ix1, iy1);

    // Bilinear blend
    double vx = (1-tx)*(1-ty)*vx00 + tx*(1-ty)*vx10
              + (1-tx)*ty   *vx01 + tx*ty    *vx11;
    double vy = (1-tx)*(1-ty)*vy00 + tx*(1-ty)*vy10
        + (1 - tx) * ty * vy01 + tx * ty * vy11;

    return { vx, vy, 0.0 };
}

} // namespace gz::waves