#include "krykov_e_sobel_op/seq/include/ops_seq.hpp"

#include <array>
#include <cmath>
#include <cstddef>
#include <vector>

#include "krykov_e_sobel_op/common/include/common.hpp"

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

  grayscale_.resize(static_cast<size_t>(width_) * static_cast<size_t>(height_));
  // RGB → grayscale
  for (int i = 0; i < width_ * height_; ++i) {
    const Pixel &p = img.data[i];
    grayscale_[i] = static_cast<int>((0.299 * p.r) + (0.587 * p.g) + (0.114 * p.b));
  }
  GetOutput().assign(static_cast<size_t>(width_) * static_cast<size_t>(height_), 0);
  return true;
}

bool KrykovESobelOpSEQ::RunImpl() {
  const std::array<std::array<int, 3>, 3> gx_kernel = {{{-1, 0, 1}, {-2, 0, 2}, {-1, 0, 1}}};

  const std::array<std::array<int, 3>, 3> gy_kernel = {{{-1, -2, -1}, {0, 0, 0}, {1, 2, 1}}};

  for (int row = 1; row < height_ - 1; ++row) {
    for (int col = 1; col < width_ - 1; ++col) {
      int gx = 0;
      int gy = 0;

      for (int ky = -1; ky <= 1; ++ky) {
        for (int kx = -1; kx <= 1; ++kx) {
          int pixel = grayscale_[((row + ky) * width_) + (col + kx)];
          gx += pixel * gx_kernel.at(ky + 1).at(kx + 1);
          gy += pixel * gy_kernel.at(ky + 1).at(kx + 1);
        }
      }

      int magnitude = static_cast<int>(std::sqrt(static_cast<double>((gx * gx) + (gy * gy))));

      GetOutput()[(row * width_) + col] = magnitude;
    }
  }

  return true;
}

bool KrykovESobelOpSEQ::PostProcessingImpl() {
  return true;
}

}  // namespace krykov_e_sobel_op
