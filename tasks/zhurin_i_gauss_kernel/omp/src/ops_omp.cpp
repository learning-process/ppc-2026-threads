#include "zhurin_i_gauss_kernel/omp/include/ops_omp.hpp"

#include <algorithm>
#include <vector>

#include "zhurin_i_gauss_kernel/common/include/common.hpp"

namespace zhurin_i_gauss_kernel {

bool ZhurinIGaussKernelOMP::ValidationImpl() {
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

bool ZhurinIGaussKernelOMP::PreProcessingImpl() {
  const auto &in = GetInput();
  width_ = std::get<0>(in);
  height_ = std::get<1>(in);
  num_parts_ = std::get<2>(in);

  padded_.assign(height_ + 2, std::vector<int>(width_ + 2, 0));
  const auto &img = std::get<3>(in);
  for (int i = 0; i < height_; ++i) {
    std::copy(img[i].begin(), img[i].end(), padded_[i + 1].begin() + 1);
  }

  result_.assign(height_, std::vector<int>(width_, 0));
  return true;
}

ZhurinIGaussKernelOMP::ZhurinIGaussKernelOMP(const InType &in) {
  SetTypeOfTask(GetStaticTypeOfTask());
  GetInput() = in;
  GetOutput() = OutType{};
}

bool ZhurinIGaussKernelOMP::RunImpl() {
  int w = width_;
  int h = height_;
  int np = num_parts_;
  int base_width = w / np;
  int remainder = w % np;

#pragma omp parallel for schedule(static) default(none) shared(w, h, base_width, remainder, np, kKernel, kShift)
  for (int part = 0; part < np; ++part) {
    int part_width = base_width + (part < remainder ? 1 : 0);
    int x_start = part * base_width + std::min(part, remainder);
    int x_end = x_start + part_width;

    for (int i = 1; i <= h; ++i) {
      for (int j = x_start + 1; j <= x_end; ++j) {
        int sum = padded_[i - 1][j - 1] * kKernel[0][0] + padded_[i - 1][j] * kKernel[0][1] +
                  padded_[i - 1][j + 1] * kKernel[0][2] + padded_[i][j - 1] * kKernel[1][0] +
                  padded_[i][j] * kKernel[1][1] + padded_[i][j + 1] * kKernel[1][2] +
                  padded_[i + 1][j - 1] * kKernel[2][0] + padded_[i + 1][j] * kKernel[2][1] +
                  padded_[i + 1][j + 1] * kKernel[2][2];
        result_[i - 1][j - 1] = sum >> kShift;
      }
    }
  }
  return true;
}

bool ZhurinIGaussKernelOMP::PostProcessingImpl() {
  GetOutput() = result_;
  return true;
}

}  // namespace zhurin_i_gauss_kernel
