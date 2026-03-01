#include "fedoseev_linear_image_filtering_vertical/seq/include/ops_seq.hpp"

#include <algorithm>
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
  if (input.data.size() != static_cast<size_t>(input.width * input.height)) {
    return false;
  }
  return true;
}

bool LinearImageFilteringVerticalSeq::PreProcessingImpl() {
  const InType &input = GetInput();
  OutType output;
  output.width = input.width;
  output.height = input.height;
  output.data.resize(input.width * input.height, 0);
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

  const int kernel[3][3] = {{1, 2, 1}, {2, 4, 2}, {1, 2, 1}};
  const int kernel_sum = 16;

  auto get_pixel = [&](int x, int y) -> int {
    x = std::clamp(x, 0, w - 1);
    y = std::clamp(y, 0, h - 1);
    return src[y * w + x];
  };

  for (int y = 0; y < h; ++y) {
    for (int x = 0; x < w; ++x) {
      int sum = 0;
      for (int ky = -1; ky <= 1; ++ky) {
        for (int kx = -1; kx <= 1; ++kx) {
          sum += get_pixel(x + kx, y + ky) * kernel[ky + 1][kx + 1];
        }
      }
      dst[y * w + x] = sum / kernel_sum;
    }
  }

  return true;
}

bool LinearImageFilteringVerticalSeq::PostProcessingImpl() {
  return true;
}

}  // namespace fedoseev_linear_image_filtering_vertical
