#include "zhurin_i_gauss_kernel_seq/seq/include/ops_seq.hpp"

#include <algorithm>
#include <cstddef>
#include <utility>
#include <vector>

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
  if (std::cmp_not_equal(img.size(), h)) {  // вместо static_cast<int>(img.size()) != h
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

  result_.assign(height_, std::vector<int>(width_, 0));

  std::vector<std::vector<int>> padded(height_ + 2, std::vector<int>(width_ + 2, 0));
  for (int i = 0; i < height_; ++i) {
    std::copy(image_[i].begin(), image_[i].end(), padded[i + 1].begin() + 1);
  }

  int basewidth_ = width_ / num_parts_;
  int remainder = width_ % num_parts_;
  int xStart = 0;

  auto convolveAt = [&](int row, int col) -> int {
    int sum = 0;
    for (int ki = 0; ki < 3; ++ki) {
      for (int kj = 0; kj < 3; ++kj) {
        sum += padded[row - 1 + ki][col - 1 + kj] * kKernel[ki][kj];
      }
    }
    return sum >> kShift;
  };

  for (int part = 0; part < num_parts_; ++part) {
    int partwidth_ = basewidth_ + (part < remainder ? 1 : 0);
    int xEnd = xStart + partwidth_;

    for (int i = 1; i <= height_; ++i) {
      for (int j = xStart + 1; j <= xEnd; ++j) {
        result_[i - 1][j - 1] = convolveAt(i, j);
      }
    }

    xStart = xEnd;
  }

  GetOutput() = std::move(result_);
  output_written_ = true;
  return true;
}

bool ZhurinIGaussKernelSEQ::RunImpl() {
  return true;
}

bool ZhurinIGaussKernelSEQ::PostProcessingImpl() {
  return output_written_;
}

}  // namespace zhurin_i_gauss_kernel_seq
