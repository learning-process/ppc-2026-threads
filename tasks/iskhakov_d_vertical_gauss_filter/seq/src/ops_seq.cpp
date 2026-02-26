#include "iskhakov_d_vertical_gauss_filter/seq/include/ops_seq.hpp"

#include <array>
#include <numeric>
#include <vector>

#include "iskhakov_d_vertical_gauss_filter/common/include/common.hpp"
#include "util/include/util.hpp"

namespace iskhakov_d_vertical_gauss_filter {

namespace {
const int DIV_CONST = 16;
const std::array<std::array<int, 3>, 3> GAUSS_KERNEL = {{{1, 2, 1}, {2, 4, 2}, {1, 2, 1}}};
}  // namespace

IskhakovDVerticalGaussFilterSEQ::IskhakovDVerticalGaussFilterSEQ(const InType &in) {
  SetTypeOfTask(GetStaticTypeOfTask());
  GetInput() = in;
  GetOutput() = OutType{};
}

bool IskhakovDVerticalGaussFilterSEQ::ValidationImpl() {
  const auto &in = GetInput();

  if (in.width <= 0 || in.height <= 0) {
    return false;
  }
  if (in.data.size() != static_cast<size_t>(in.width * in.height)) {
    return false;
  }
  return true;
}

bool IskhakovDVerticalGaussFilterSEQ::PreProcessingImpl() {
  return true;
}

static uint8_t IskhakovDGetPixelMirrorSeq(const std::vector<uint8_t> &res, int col, int row, int width, int height) {
  if (col < 0) {
    col = -col - 1;
  } else if (col >= width) {
    col = 2 * width - col - 1;
  }
  if (row < 0) {
    row = -row - 1;
  } else if (row >= height) {
    row = 2 * height - row - 1;
  }
  return res[row * width + col];
}

bool IskhakovDVerticalGaussFilterSEQ::RunImpl() {
  const auto &in = GetInput();

  int width = in.width;
  int height = in.height;
  const std::vector<uint8_t> &matrix = in.data;
  std::vector<uint8_t> result(width * height);

  for (int horizontal = 0; horizontal < width; ++horizontal) {
    for (int vertical = 0; vertical < height; ++vertical) {
      int sum = 0;

      sum += 1 * IskhakovDGetPixelMirrorSeq(matrix, horizontal - 1, vertical - 1, width, height);
      sum += 2 * IskhakovDGetPixelMirrorSeq(matrix, horizontal, vertical - 1, width, height);
      sum += 1 * IskhakovDGetPixelMirrorSeq(matrix, horizontal + 1, vertical - 1, width, height);

      sum += 2 * IskhakovDGetPixelMirrorSeq(matrix, horizontal - 1, vertical, width, height);
      sum += 4 * IskhakovDGetPixelMirrorSeq(matrix, horizontal, vertical, width, height);
      sum += 2 * IskhakovDGetPixelMirrorSeq(matrix, horizontal + 1, vertical, width, height);

      sum += 1 * IskhakovDGetPixelMirrorSeq(matrix, horizontal - 1, vertical + 1, width, height);
      sum += 2 * IskhakovDGetPixelMirrorSeq(matrix, horizontal, vertical + 1, width, height);
      sum += 1 * IskhakovDGetPixelMirrorSeq(matrix, horizontal + 1, vertical + 1, width, height);

      result[vertical * width + horizontal] = static_cast<uint8_t>(sum / DIV_CONST);
    }
  }
  GetOutput().width = width;
  GetOutput().height = height;
  GetOutput().data = std::move(result);
  return true;
}

bool IskhakovDVerticalGaussFilterSEQ::PostProcessingImpl() {
  return true;
}

}  // namespace iskhakov_d_vertical_gauss_filter
