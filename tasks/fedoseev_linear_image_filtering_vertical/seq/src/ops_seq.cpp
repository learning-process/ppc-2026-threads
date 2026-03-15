#include "fedoseev_linear_image_filtering_vertical/seq/include/ops_seq.hpp"

#include <algorithm>
#include <array>
#include <cstddef>
#include <vector>

#include "fedoseev_linear_image_filtering_vertical/common/include/common.hpp"

namespace fedoseev_linear_image_filtering_vertical {

LinearImageFilteringVerticalSeq::LinearImageFilteringVerticalSeq(const InType &in) {
  SetTypeOfTask(GetStaticTypeOfTask());
  GetInput() = in;
  GetOutput() = InType{};
}

bool LinearImageFilteringVerticalSeq::ValidationImpl() {
  const InType &input = GetInput();
  if (input.width < 3 || input.height < 3) {
    return false;
  }
  if (input.data.size() != static_cast<size_t>(input.width) * static_cast<size_t>(input.height)) {
    return false;
  }
  return true;
}

bool LinearImageFilteringVerticalSeq::PreProcessingImpl() {
  const InType &input = GetInput();
  OutType output;
  output.width = input.width;
  output.height = input.height;
  output.data.resize(static_cast<size_t>(input.width) * static_cast<size_t>(input.height), 0);
  GetOutput() = output;
  return true;
}

bool LinearImageFilteringVerticalSeq::RunImpl() {
  const InType &input = GetInput();
  OutType &output = GetOutput();

  int w = input.width;
  int h = input.height;
  const std::vector<int> &src = input.data;
  std::vector<int> &dst = output.data;

  const std::array<std::array<int, 3>, 3> kernel = {{{{1, 2, 1}}, {{2, 4, 2}}, {{1, 2, 1}}}};
  const int kernel_sum = 16;

  auto get_pixel = [&](int col, int row) -> int {
    col = std::clamp(col, 0, w - 1);
    row = std::clamp(row, 0, h - 1);
    return src[(static_cast<size_t>(row) * static_cast<size_t>(w)) + static_cast<size_t>(col)];
  };

  for (int row = 0; row < h; ++row) {
    for (int col = 0; col < w; ++col) {
      int sum = 0;
      for (int ky = -1; ky <= 1; ++ky) {
        for (int kx = -1; kx <= 1; ++kx) {
          sum += get_pixel(col + kx, row + ky) * kernel.at(ky + 1).at(kx + 1);
        }
      }
      dst[(static_cast<size_t>(row) * static_cast<size_t>(w)) + static_cast<size_t>(col)] = sum / kernel_sum;
    }
  }

  return true;
}

bool LinearImageFilteringVerticalSeq::PostProcessingImpl() {
  return true;
}

}  // namespace fedoseev_linear_image_filtering_vertical
