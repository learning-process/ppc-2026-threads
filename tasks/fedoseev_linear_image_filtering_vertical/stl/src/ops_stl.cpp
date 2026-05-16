#include "fedoseev_linear_image_filtering_vertical/stl/include/ops_stl.hpp"

#include <algorithm>
#include <array>
#include <cstddef>
#include <thread>
#include <vector>

#include "fedoseev_linear_image_filtering_vertical/common/include/common.hpp"

namespace fedoseev_linear_image_filtering_vertical {

namespace {
int GetPixel(const std::vector<int> &src, int w, int h, int col, int row) {
  col = std::clamp(col, 0, w - 1);
  row = std::clamp(row, 0, h - 1);
  return src[(static_cast<size_t>(row) * static_cast<size_t>(w)) + static_cast<size_t>(col)];
}

int ComputePixel(const std::vector<int> &src, int w, int h, int col, int row,
                 const std::array<std::array<int, 3>, 3> &kernel, int kernel_sum) {
  int sum = 0;
  for (int ky = 0; ky < 3; ++ky) {
    for (int kx = 0; kx < 3; ++kx) {
      int px = col + kx - 1;
      int py = row + ky - 1;
      sum += GetPixel(src, w, h, px, py) * kernel.at(ky).at(kx);
    }
  }
  return sum / kernel_sum;
}

void ProcessBlock(const std::vector<int> &src, std::vector<int> &dst, int w, int h, int col_start, int col_end,
                  const std::array<std::array<int, 3>, 3> &kernel, int kernel_sum) {
  for (int row = 0; row < h; ++row) {
    for (int col = col_start; col < col_end; ++col) {
      dst[(static_cast<size_t>(row) * static_cast<size_t>(w)) + static_cast<size_t>(col)] =
          ComputePixel(src, w, h, col, row, kernel, kernel_sum);
    }
  }
}
}  // namespace

LinearImageFilteringVerticalSTL::LinearImageFilteringVerticalSTL(const InType &in) {
  SetTypeOfTask(GetStaticTypeOfTask());
  GetInput() = in;
  GetOutput() = InType{};
}

bool LinearImageFilteringVerticalSTL::ValidationImpl() {
  const InType &input = GetInput();
  if (input.width < 3 || input.height < 3) {
    return false;
  }
  if (input.data.size() != static_cast<size_t>(input.width) * static_cast<size_t>(input.height)) {
    return false;
  }
  return true;
}

bool LinearImageFilteringVerticalSTL::PreProcessingImpl() {
  const InType &input = GetInput();
  OutType output;
  output.width = input.width;
  output.height = input.height;
  output.data.resize(static_cast<size_t>(input.width) * static_cast<size_t>(input.height), 0);
  GetOutput() = output;
  return true;
}

bool LinearImageFilteringVerticalSTL::RunImpl() {
  const InType &input = GetInput();
  OutType &output = GetOutput();

  int w = input.width;
  int h = input.height;
  const std::vector<int> &src = input.data;
  std::vector<int> &dst = output.data;

  const std::array<std::array<int, 3>, 3> kernel = {{{{1, 2, 1}}, {{2, 4, 2}}, {{1, 2, 1}}}};
  const int block_width = 64;
  const int num_blocks = (w + block_width - 1) / block_width;

  std::vector<std::thread> threads;
  threads.reserve(static_cast<size_t>(num_blocks));

  for (int block_idx = 0; block_idx < num_blocks; ++block_idx) {
    int col_start = block_idx * block_width;
    int col_end = std::min(col_start + block_width, w);
    threads.emplace_back([&src, &dst, w, h, col_start, col_end, &kernel]() {
      ProcessBlock(src, dst, w, h, col_start, col_end, kernel, 16);
    });
  }

  for (auto &t : threads) {
    t.join();
  }
  return true;
}

bool LinearImageFilteringVerticalSTL::PostProcessingImpl() {
  return true;
}

}  // namespace fedoseev_linear_image_filtering_vertical
