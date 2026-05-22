#include "belov_e_sobel/omp/include/ops_omp.hpp"

#include <omp.h>

#include <algorithm>
#include <cmath>
#include <vector>

#include "belov_e_sobel/common/include/common.hpp"
#include "util/include/util.hpp"

namespace belov_e_sobel {

BelovESobelOMP::BelovESobelOMP(const InType &in) {
  SetTypeOfTask(GetStaticTypeOfTask());
  GetInput() = in;
  GetOutput() = in;
}

bool BelovESobelOMP::ValidationImpl() {
  return !std::get<0>(GetInput()).empty() && (std::get<1>(GetInput()) > 0) && (std::get<2>(GetInput()) > 0);
}

bool BelovESobelOMP::PreProcessingImpl() {
  return true;
}

bool BelovESobelOMP::RunImpl() {
  const std::vector<uint8_t> &input = std::get<0>(GetInput());
  std::vector<uint8_t> &output = std::get<0>(GetOutput());
  int width = std::get<1>(GetInput());
  int height = std::get<2>(GetInput());

#pragma omp parallel for schedule(static)
  for (int y = 0; y < height; ++y) {
    for (int x = 0; x < width; ++x) {
      int x_minus = std::clamp(x - 1, 0, width - 1);
      int x_plus = std::clamp(x + 1, 0, width - 1);
      int y_minus = std::clamp(y - 1, 0, height - 1);
      int y_plus = std::clamp(y + 1, 0, height - 1);

      float gx = (-1.0f * input[y_minus * width + x_minus]) + (1.0f * input[y_minus * width + x_plus]) +
                 (-2.0f * input[y * width + x_minus]) + (2.0f * input[y * width + x_plus]) +
                 (-1.0f * input[y_plus * width + x_minus]) + (1.0f * input[y_plus * width + x_plus]);

      float gy = (-1.0f * input[y_minus * width + x_minus]) - (2.0f * input[y_minus * width + x]) -
                 (1.0f * input[y_minus * width + x_plus]) + (1.0f * input[y_plus * width + x_minus]) +
                 (2.0f * input[y_plus * width + x]) + (1.0f * input[y_plus * width + x_plus]);

      float magnitude = std::sqrt(gx * gx + gy * gy);
      output[y * width + x] = static_cast<uint8_t>(std::min(255.0f, magnitude));
    }
  }

  return true;
}

bool BelovESobelOMP::PostProcessingImpl() {
  return !std::get<0>(GetOutput()).empty() && (std::get<1>(GetOutput()) > 0) && (std::get<2>(GetOutput()) > 0);
}

}  // namespace belov_e_sobel
