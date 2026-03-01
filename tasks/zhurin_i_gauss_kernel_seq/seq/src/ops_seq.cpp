#include "zhurin_i_gauss_kernel_seq/seq/include/ops_seq.hpp"

#include <algorithm>
#include <utility>
#include <vector>

#include "zhurin_i_gauss_kernel_seq/common/include/common.hpp"

namespace zhurin_i_gauss_kernel_seq {

ZhurinIGaussKernelSEQ::ZhurinIGaussKernelSEQ(const InType &in) {
  SetTypeOfTask(GetStaticTypeOfTask());
  GetInput() = in;
  GetOutput() = OutType{};
}

bool ZhurinIGaussKernelSEQ::ValidationImpl() {
  const auto &in = GetInput();
  int w = std::get<0>(in);
  int h = std::get<1>(in);
  int parts = std::get<2>(in);
  const auto &img = std::get<3>(in);

  if (w <= 0 || h <= 0 || parts <= 0 || parts > w) {
    return false;
  }
  if (std::cmp_not_equal(img.size(), h)) {
    return false;
  }
  for (int i = 0; i < h; ++i) {
    if (std::cmp_not_equal(img[i].size(), w)) {
      return false;
    }
  }
  return true;
}

bool ZhurinIGaussKernelSEQ::PreProcessingImpl() {
  const auto &in = GetInput();
  width_ = std::get<0>(in);
  height_ = std::get<1>(in);
  num_parts_ = std::get<2>(in);
  image_ = std::get<3>(in);

  padded_.assign(height_ + 2, std::vector<int>(width_ + 2, 0));
  for (int i = 0; i < height_; ++i) {
    std::copy(image_[i].begin(), image_[i].end(), padded_[i + 1].begin() + 1);
  }

  result_.assign(height_, std::vector<int>(width_, 0));
  output_written_ = false;
  return true;
}

bool ZhurinIGaussKernelSEQ::RunImpl() {
  int base_width = width_ / num_parts_;
  int remainder = width_ % num_parts_;
  int x_start = 0;

  auto convolve_at = [&](int row, int col) -> int {
    int sum = 0;
    sum += padded_[row - 1][col - 1] * kKernel[0][0];
    sum += padded_[row - 1][col] * kKernel[0][1];
    sum += padded_[row - 1][col + 1] * kKernel[0][2];
    sum += padded_[row][col - 1] * kKernel[1][0];
    sum += padded_[row][col] * kKernel[1][1];
    sum += padded_[row][col + 1] * kKernel[1][2];
    sum += padded_[row + 1][col - 1] * kKernel[2][0];
    sum += padded_[row + 1][col] * kKernel[2][1];
    sum += padded_[row + 1][col + 1] * kKernel[2][2];
    return sum >> kShift;
  };

  for (int part = 0; part < num_parts_; ++part) {
    int part_width = base_width + (part < remainder ? 1 : 0);
    int x_end = x_start + part_width;

    for (int i = 1; i <= height_; ++i) {
      for (int j = x_start + 1; j <= x_end; ++j) {
        result_[i - 1][j - 1] = convolve_at(i, j);
      }
    }
    x_start = x_end;
  }
  return true;
}

bool ZhurinIGaussKernelSEQ::PostProcessingImpl() {
  GetOutput() = result_; 
  output_written_ = true;
  return true;
}

}  // namespace zhurin_i_gauss_kernel_seq
