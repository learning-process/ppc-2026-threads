#include "fedoseev_linear_image_filtering_vertical/omp/include/ops_omp.hpp"

#include <omp.h>

#include <algorithm>
#include <array>
#include <cstddef>
#include <vector>

#include "fedoseev_linear_image_filtering_vertical/common/include/common.hpp"

namespace fedoseev_linear_image_filtering_vertical {

namespace {
int GetPixel(const std::vector<int> &src, int w, int h, int col, int row) {
  col = std::clamp(col, 0, w - 1);
  row = std::clamp(row, 0, h - 1);
  return src[(static_cast<size_t>(row) * static_cast<size_t>(w)) + static_cast<size_t>(col)];
}
}  // namespace

LinearImageFilteringVerticalOMP::LinearImageFilteringVerticalOMP(const InType &in) {
  SetTypeOfTask(GetStaticTypeOfTask());
  GetInput() = in;
  GetOutput() = InType{};
}

bool LinearImageFilteringVerticalOMP::ValidationImpl() {
  const InType &input = GetInput();
  if (input.width < 3 || input.height < 3) {
    return false;
  }
  if (input.data.size() != static_cast<size_t>(input.width) * static_cast<size_t>(input.height)) {
    return false;
  }
  return true;
}

bool LinearImageFilteringVerticalOMP::PreProcessingImpl() {
  const InType &input = GetInput();
  OutType output;
  output.width = input.width;
  output.height = input.height;
  output.data.resize(static_cast<size_t>(input.width) * static_cast<size_t>(input.height), 0);
  GetOutput() = output;
  return true;
}

bool LinearImageFilteringVerticalOMP::RunImpl() {
  const InType &input = GetInput();
  OutType &output = GetOutput();

  int w = input.width;
  int h = input.height;
  const std::vector<int> &src = input.data;
  std::vector<int> &dst = output.data;

  const std::array<std::array<int, 3>, 3> kernel = {{{{1, 2, 1}}, {{2, 4, 2}}, {{1, 2, 1}}}};
  const int kernel_sum = 16;
  const int block_width = 64;

#pragma omp parallel for schedule(static) default(none) shared(w, h, src, dst, kernel, kernel_sum, block_width)
  for (int col_start = 0; col_start < w; col_start += block_width) {
    int col_end = std::min(col_start + block_width, w);
    for (int row = 0; row < h; ++row) {
      for (int col = col_start; col < col_end; ++col) {
        int sum = 0;
        for (int ky = 0; ky < 3; ++ky) {
          for (int kx = 0; kx < 3; ++kx) {
            int px = col + kx - 1;
            int py = row + ky - 1;
            int pixel = GetPixel(src, w, h, px, py);
            sum += pixel * kernel.at(ky).at(kx);
          }
        }
        dst[(static_cast<size_t>(row) * static_cast<size_t>(w)) + static_cast<size_t>(col)] = sum / kernel_sum;
      }
    }
  }

  return true;
}

bool LinearImageFilteringVerticalOMP::PostProcessingImpl() {
  return true;
}

}  // namespace fedoseev_linear_image_filtering_vertical
