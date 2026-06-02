#include "shakirova_e_sobel_edge_detection/omp/include/ops_omp.hpp"

#include <algorithm>
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
  auto &out = GetOutput();
  const int h = height_;
  const int w = width_;
  const int *inp = input_.data();

#pragma omp parallel for default(none) shared(out, inp) firstprivate(h, w) schedule(static)
  for (int row = 1; row < h - 1; ++row) {
    const int *prev = inp + static_cast<ptrdiff_t>(row - 1) * w;
    const int *curr = inp + static_cast<ptrdiff_t>(row) * w;
    const int *next = inp + static_cast<ptrdiff_t>(row + 1) * w;

    for (int col = 1; col < w - 1; ++col) {
      const int gx =
          -prev[col - 1] + prev[col + 1] - (2 * curr[col - 1]) + (2 * curr[col + 1]) - next[col - 1] + next[col + 1];

      const int gy = -prev[col - 1] - (2 * prev[col]) - prev[col + 1] + next[col - 1] + (2 * next[col]) + next[col + 1];

      const int abs_gx = std::abs(gx);
      const int abs_gy = std::abs(gy);
      const int magnitude = (std::max(abs_gx, abs_gy) * 123 + std::min(abs_gx, abs_gy) * 51) >> 7;
      out[(row * w) + col] = std::min(magnitude, 255);
    }
  }

  return true;
}

bool ShakirovaESobelEdgeDetectionOMP::PostProcessingImpl() {
  return true;
}

}  // namespace shakirova_e_sobel_edge_detection
