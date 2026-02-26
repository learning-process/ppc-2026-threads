#include "shakirova_e_sobel_edge_detection/seq/include/ops_seq.hpp"

#include <numeric>
#include <vector>

#include "shakirova_e_sobel_edge_detection/common/include/common.hpp"
#include "util/include/util.hpp"

namespace shakirova_e_sobel_edge_detection {

ShakirovaESobelEdgeDetectionSEQ::ShakirovaESobelEdgeDetectionSEQ(const InType &in) {
  SetTypeOfTask(GetStaticTypeOfTask());
  GetInput() = in;
  GetOutput().clear();
}

bool ShakirovaESobelEdgeDetectionSEQ::ValidationImpl() {
  return GetInput().IsValid();
}

bool ShakirovaESobelEdgeDetectionSEQ::PreProcessingImpl() {
  const auto &img = GetInput();
  width_  = img.width;
  height_ = img.height;
  input_  = img.pixels;

  GetOutput().assign(static_cast<size_t>(width_ * height_), 0);
  return true;
}

bool ShakirovaESobelEdgeDetectionSEQ::RunImpl() {
  const std::array<std::array<int, 3>, 3> kGx = {{{-1, 0, 1}, {-2, 0, 2}, {-1, 0, 1}}};
  const std::array<std::array<int, 3>, 3> kGy = {{{-1, -2, -1}, {0, 0, 0}, {1, 2, 1}}};

  auto &out = GetOutput();

  for (int y = 1; y < height_ - 1; ++y) {
    for (int x = 1; x < width_ - 1; ++x) {
      int gx = 0;
      int gy = 0;

      for (int ky = -1; ky <= 1; ++ky) {
        for (int kx = -1; kx <= 1; ++kx) {
          const int pixel = input_[((y + ky) * width_) + (x + kx)];
          gx += pixel * kGx[ky + 1][kx + 1];
          gy += pixel * kGy[ky + 1][kx + 1];
        }
      }

      const int magnitude =
          static_cast<int>(std::sqrt(static_cast<double>(gx * gx + gy * gy)));
      out[(y * width_) + x] = std::clamp(magnitude, 0, 255);
    }
  }

  return true;
}

bool ShakirovaESobelEdgeDetectionSEQ::PostProcessingImpl() {
  return true;
}

}  // namespace shakirova_e_sobel_edge_detection
