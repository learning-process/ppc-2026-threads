#include "krykov_e_sobel_op/seq/include/ops_seq.hpp"

#include <numeric>
#include <vector>

#include "krykov_e_sobel_op/common/include/common.hpp"
#include "util/include/util.hpp"

namespace krykov_e_sobel_op {

KrykovESobelOpSEQ::KrykovESobelOpSEQ(const InType &in) {
  SetTypeOfTask(GetStaticTypeOfTask());
  GetInput() = in;
  GetOutput().clear();
}

bool KrykovESobelOpSEQ::ValidationImpl() {
  const auto &img = GetInput();
  return img.width > 2 && img.height > 2 && static_cast<int>(img.data.size()) == img.width * img.height;
}

bool KrykovESobelOpSEQ::PreProcessingImpl() {
  const auto &img = GetInput();

  width_ = img.width;
  height_ = img.height;

  grayscale_.resize(width_ * height_);
  // RGB → grayscale
  for (int i = 0; i < width_ * height_; ++i) {
    const Pixel &p = img.data[i];
    grayscale_[i] = static_cast<int>(0.299 * p.r + 0.587 * p.g + 0.114 * p.b);
  }
  GetOutput().assign(width_ * height_, 0);
  return true;
}

bool KrykovESobelOpSEQ::RunImpl() {
  const int gx_kernel[3][3] = {{-1, 0, 1}, {-2, 0, 2}, {-1, 0, 1}};

  const int gy_kernel[3][3] = {{-1, -2, -1}, {0, 0, 0}, {1, 2, 1}};

  for (int y = 1; y < height_ - 1; ++y) {
    for (int x = 1; x < width_ - 1; ++x) {
      int gx = 0;
      int gy = 0;

      for (int ky = -1; ky <= 1; ++ky) {
        for (int kx = -1; kx <= 1; ++kx) {
          int pixel = grayscale_[(y + ky) * width_ + (x + kx)];
          gx += pixel * gx_kernel[ky + 1][kx + 1];
          gy += pixel * gy_kernel[ky + 1][kx + 1];
        }
      }

      int magnitude = static_cast<int>(std::sqrt(static_cast<double>(gx * gx + gy * gy)));

      GetOutput()[y * width_ + x] = magnitude;
    }
  }

  return true;
}

bool KrykovESobelOpSEQ::PostProcessingImpl() {
  return true;
}

}  // namespace krykov_e_sobel_op
