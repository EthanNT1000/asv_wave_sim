// gz-waves/include/gz/waves/WaterCurrentGrid.hh
#pragma once

#include <string>
#include <vector>
#include "gz/waves/CGALTypes.hh"

namespace gz::waves {

/// Loads and samples a pre-processed binary water current grid
/// produced by preprocess_hecras.py (format WCRG v1).
class WaterCurrentGrid
{
public:
    WaterCurrentGrid() = default;

    /// Load grid from binary file. Returns false and logs error on failure.
    bool LoadFromFile(const std::string& bin_path);

    /// True if a grid has been successfully loaded.
    bool IsLoaded() const { return loaded_; }

    /// Bilinear interpolation of water current velocity at world position (x, y).
    /// Returns (0,0,0) if (x,y) is outside the grid extent.
    gz::cgal::Vector3 SampleAt(double x, double y) const;

    /// Grid extents in local frame (for debugging / visualisation).
    double XMin() const { return x_min_; }
    double YMin() const { return y_min_; }
    double XMax() const { return x_min_ + cell_size_x_ * (nx_ - 1); }
    double YMax() const { return y_min_ + cell_size_y_ * (ny_ - 1); }

private:
    bool loaded_      = false;
    int  nx_          = 0;
    int  ny_          = 0;
    double x_min_     = 0.0;
    double y_min_     = 0.0;
    double cell_size_x_ = 1.0;
    double cell_size_y_ = 1.0;

    // Interleaved (vx, vy) pairs, row-major: index = iy * nx_ + ix
    std::vector<float> data_;   // size = nx_ * ny_ * 2
};

} // namespace gz::waves
