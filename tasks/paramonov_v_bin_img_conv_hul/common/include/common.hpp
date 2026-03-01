#pragma once

#include <cstdint>
#include <vector>

#include "task/include/task.hpp"

namespace paramonov_v_bin_img_conv_hul {

struct PixelPoint {
  int row{0};
  int col{0};

  PixelPoint() = default;
  PixelPoint(int r, int c) : row(r), col(c) {}

  bool operator==(const PixelPoint &other) const {
    return row == other.row && col == other.col;
  }

  bool operator<(const PixelPoint &other) const {
    return (row == other.row) ? col < other.col : row < other.row;
  }
};

struct GrayImage {
  std::vector<uint8_t> pixels;
  int rows = 0;
  int cols = 0;
};

using InputType = GrayImage;
using OutputType = std::vector<std::vector<PixelPoint>>;
using HullTaskBase = ppc::task::Task<InputType, OutputType>;

}  // namespace paramonov_v_bin_img_conv_hul
