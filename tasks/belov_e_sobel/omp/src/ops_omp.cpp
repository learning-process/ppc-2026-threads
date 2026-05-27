#include "belov_e_sobel/omp/include/ops_omp.hpp"

#include <omp.h>

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <vector>

#include "belov_e_sobel/common/include/common.hpp"

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

#pragma omp parallel for default(none) shared(input, output, width, height) schedule(static)
  for (int row = 0; row < height; ++row) {
    for (int col = 0; col < width; ++col) {
      int x_minus = std::clamp(col - 1, 0, width - 1);
      int x_plus = std::clamp(col + 1, 0, width - 1);
      int y_minus = std::clamp(row - 1, 0, height - 1);
      int y_plus = std::clamp(row + 1, 0, height - 1);

      float gx = (-1.0F * static_cast<float>(input[(static_cast<std::size_t>(y_minus) * width) + x_minus])) +
                 (1.0F * static_cast<float>(input[(static_cast<std::size_t>(y_minus) * width) + x_plus])) +
                 (-2.0F * static_cast<float>(input[(static_cast<std::size_t>(row) * width) + x_minus])) +
                 (2.0F * static_cast<float>(input[(static_cast<std::size_t>(row) * width) + x_plus])) +
                 (-1.0F * static_cast<float>(input[(static_cast<std::size_t>(y_plus) * width) + x_minus])) +
                 (1.0F * static_cast<float>(input[(static_cast<std::size_t>(y_plus) * width) + x_plus]));

      float gy = (-1.0F * static_cast<float>(input[(static_cast<std::size_t>(y_minus) * width) + x_minus])) -
                 (2.0F * static_cast<float>(input[(static_cast<std::size_t>(y_minus) * width) + col])) -
                 (1.0F * static_cast<float>(input[(static_cast<std::size_t>(y_minus) * width) + x_plus])) +
                 (1.0F * static_cast<float>(input[(static_cast<std::size_t>(y_plus) * width) + x_minus])) +
                 (2.0F * static_cast<float>(input[(static_cast<std::size_t>(y_plus) * width) + col])) +
                 (1.0F * static_cast<float>(input[(static_cast<std::size_t>(y_plus) * width) + x_plus]));

      float magnitude = std::sqrt((gx * gx) + (gy * gy));
      output[(static_cast<std::size_t>(row) * width) + col] = static_cast<uint8_t>(std::min(255.0F, magnitude));
    }
  }

  return true;
}

bool BelovESobelOMP::PostProcessingImpl() {
  return !std::get<0>(GetOutput()).empty() && (std::get<1>(GetOutput()) > 0) && (std::get<2>(GetOutput()) > 0);
}

}  // namespace belov_e_sobel
