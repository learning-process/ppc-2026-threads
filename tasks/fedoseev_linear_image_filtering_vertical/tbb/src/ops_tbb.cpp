#include "fedoseev_linear_image_filtering_vertical/tbb/include/ops_tbb.hpp"

#include <tbb/blocked_range.h>
#include <tbb/parallel_for.h>

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

LinearImageFilteringVerticalTBB::LinearImageFilteringVerticalTBB(const InType &in) {
  SetTypeOfTask(GetStaticTypeOfTask());
  GetInput() = in;
  GetOutput() = InType{};
}

bool LinearImageFilteringVerticalTBB::ValidationImpl() {
  const InType &input = GetInput();
  if (input.width < 3 || input.height < 3) {
    return false;
  }
  return input.data.size() == static_cast<size_t>(input.width) * static_cast<size_t>(input.height);
}

bool LinearImageFilteringVerticalTBB::PreProcessingImpl() {
  const InType &input = GetInput();
  OutType output;
  output.width = input.width;
  output.height = input.height;
  output.data.resize(static_cast<size_t>(input.width) * static_cast<size_t>(input.height), 0);
  GetOutput() = output;
  return true;
}

bool LinearImageFilteringVerticalTBB::RunImpl() {
  const InType &input = GetInput();
  OutType &output = GetOutput();

  int w = input.width;
  int h = input.height;
  const std::vector<int> &src = input.data;
  std::vector<int> &dst = output.data;

  const std::array<std::array<int, 3>, 3> kernel = {{{{1, 2, 1}}, {{2, 4, 2}}, {{1, 2, 1}}}};
  const int kernel_sum = 16;
  const int block_width = 64;

  tbb::parallel_for(tbb::blocked_range<int>(0, w, block_width), [&](const tbb::blocked_range<int> &range) {
    for (int row = 0; row < h; ++row) {
      for (int col = range.begin(); col != range.end(); ++col) {
        int sum = 0;
        for (int ky = 0; ky < 3; ++ky) {
          for (int kx = 0; kx < 3; ++kx) {
            int px = col + kx - 1;
            int py = row + ky - 1;
            sum += GetPixel(src, w, h, px, py) * kernel.at(ky).at(kx);
          }
        }
        dst[(static_cast<size_t>(row) * static_cast<size_t>(w)) + static_cast<size_t>(col)] = sum / kernel_sum;
      }
    }
  });

  return true;
}

bool LinearImageFilteringVerticalTBB::PostProcessingImpl() {
  return true;
}

}  // namespace fedoseev_linear_image_filtering_vertical
