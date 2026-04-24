#pragma once

#include <cstdint>
#include <vector>

#include "task/include/task.hpp"

namespace dorogin_v_bin_img_conv_hull_tbb {

struct PixelPoint {
  int row{0};
  int col{0};

  PixelPoint() = default;
  PixelPoint(int row, int col) : row(row), col(col) {}

  bool operator==(const PixelPoint &other) const {
    return row == other.row && col == other.col;
  }
};

struct GrayImage {
  std::vector<uint8_t> pixels;
  int rows{0};
  int cols{0};
};

using InputType = GrayImage;
using OutputType = std::vector<std::vector<PixelPoint>>;
using BaseTask = ppc::task::Task<InputType, OutputType>;

}  // namespace dorogin_v_bin_img_conv_hull_tbb
