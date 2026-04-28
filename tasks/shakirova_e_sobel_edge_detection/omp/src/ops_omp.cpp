#include "shakirova_e_sobel_edge_detection/omp/include/ops_omp.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <vector>

#include "shakirova_e_sobel_edge_detection/common/include/common.hpp"

namespace shakirova_e_sobel_edge_detection {

ShakirovaESobelEdgeDetectionOMP::ShakirovaESobelEdgeDetectionOMP(const InType &in) {
  SetTypeOfTask(GetStaticTypeOfTask());
  GetInput() = in;
  GetOutput().clear();
}

bool ShakirovaESobelEdgeDetectionOMP::ValidationImpl() {
  return GetInput().IsValid();
}

bool ShakirovaESobelEdgeDetectionOMP::PreProcessingImpl() {
  const auto &img = GetInput();
  width_ = img.width;
  height_ = img.height;
  input_ = img.pixels;

  GetOutput().assign(static_cast<size_t>(width_) * static_cast<size_t>(height_), 0);
  return true;
}

bool ShakirovaESobelEdgeDetectionOMP::RunImpl() {
  const std::array<std::array<int, 3>, 3> k_gx = {{{-1, 0, 1}, {-2, 0, 2}, {-1, 0, 1}}};
  const std::array<std::array<int, 3>, 3> k_gy = {{{-1, -2, -1}, {0, 0, 0}, {1, 2, 1}}};

  auto &out = GetOutput();
  const int h = height_;
  const int w = width_;

#pragma omp parallel for default(none) shared(out, k_gx, k_gy) firstprivate(h, w) schedule(static)
  for (int row = 1; row < h - 1; ++row) {
    for (int col = 1; col < w - 1; ++col) {
      int gx = 0;
      int gy = 0;

      for (int ky = -1; ky <= 1; ++ky) {
        for (int kx = -1; kx <= 1; ++kx) {
          const int pixel = input_[((row + ky) * w) + (col + kx)];
          const auto ky_idx = static_cast<size_t>(ky) + 1;
          const auto kx_idx = static_cast<size_t>(kx) + 1;
          gx += pixel * k_gx.at(ky_idx).at(kx_idx);
          gy += pixel * k_gy.at(ky_idx).at(kx_idx);
        }
      }

      const int magnitude = static_cast<int>(std::sqrt(static_cast<double>((gx * gx) + (gy * gy))));
      out[(row * w) + col] = std::clamp(magnitude, 0, 255);
    }
  }

  return true;
}

bool ShakirovaESobelEdgeDetectionOMP::PostProcessingImpl() {
  return true;
}

}  // namespace shakirova_e_sobel_edge_detection
