#include "zhurin_i_gauss_kernel/omp/include/ops_omp.hpp"

#include <algorithm>
#include <utility>
#include <vector>

#include "zhurin_i_gauss_kernel/common/include/common.hpp"

#ifdef _OPENMP
#  include <omp.h>
#endif

namespace zhurin_i_gauss_kernel {

ZhurinIGaussKernelOMP::ZhurinIGaussKernelOMP(const InType &in) {
  SetTypeOfTask(GetStaticTypeOfTask());
  GetInput() = in;
  GetOutput() = OutType{};
}

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
  image_ = std::get<3>(in);

  padded_.assign(height_ + 2, std::vector<int>(width_ + 2, 0));
  for (int i = 0; i < height_; ++i) {
    std::copy(image_[i].begin(), image_[i].end(), padded_[i + 1].begin() + 1);
  }

  result_.assign(height_, std::vector<int>(width_, 0));
  output_written_ = false;
  return true;
}

bool ZhurinIGaussKernelOMP::RunImpl() {
#pragma omp parallel for schedule(static)
  for (int i = 0; i < height_; ++i) {
    for (int j = 0; j < width_; ++j) {
      int sum = 0;
      for (int ki = 0; ki < 3; ++ki) {
        for (int kj = 0; kj < 3; ++kj) {
          sum += padded_[i + ki][j + kj] * kKernel[ki][kj];
        }
      }
      result_[i][j] = sum >> kShift;
    }
  }
  return true;
}

bool ZhurinIGaussKernelOMP::PostProcessingImpl() {
  GetOutput() = std::move(result_);
  output_written_ = true;
  return true;
}

}  // namespace zhurin_i_gauss_kernel
