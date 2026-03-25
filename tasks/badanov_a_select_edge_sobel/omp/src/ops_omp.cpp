#include "badanov_a_select_edge_sobel/omp/include/ops_omp.hpp"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <vector>

#include "badanov_a_select_edge_sobel/common/include/common.hpp"

#ifdef _OPENMP
#  include <omp.h>
#endif

namespace badanov_a_select_edge_sobel {

BadanovASelectEdgeSobelOMP::BadanovASelectEdgeSobelOMP(const InType &in) {
  SetTypeOfTask(GetStaticTypeOfTask());
  GetInput() = in;
  GetOutput() = std::vector<uint8_t>();
}

bool BadanovASelectEdgeSobelOMP::ValidationImpl() {
  const auto &input = GetInput();
  return !input.empty();
}

bool BadanovASelectEdgeSobelOMP::PreProcessingImpl() {
  const auto &input = GetInput();

  width_ = static_cast<int>(std::sqrt(input.size()));
  height_ = width_;

  if (width_ * height_ != static_cast<int>(input.size())) {
    width_ = static_cast<int>(input.size());
    height_ = 1;
  }

  GetOutput() = std::vector<uint8_t>(input.size(), 0);

  return true;
}

void BadanovASelectEdgeSobelOMP::ApplySobelOperator(const std::vector<uint8_t> &input, std::vector<float> &magnitude,
                                                    float &max_magnitude) {
  max_magnitude = 0.0F;
  const int height = height_;
  const int width = width_;

  const auto &kernels_x = kKernelX;
  const auto &kernels_y = kKernelY;

#pragma omp parallel default(none) shared(input, magnitude, max_magnitude, height, width, kernels_x, kernels_y)
  {
    float local_max_magnitude = 0.0F;

#pragma omp for schedule(static)
    for (int row = 1; row < height - 1; ++row) {
      for (int col = 1; col < width - 1; ++col) {
        float gradient_x = 0.0F;
        float gradient_y = 0.0F;

        ComputeGradientAtPixel(input, row, col, gradient_x, gradient_y);

        const float magnitude_value = std::sqrt((gradient_x * gradient_x) + (gradient_y * gradient_y));
        const size_t idx = (static_cast<size_t>(row) * static_cast<size_t>(width)) + static_cast<size_t>(col);
        magnitude[idx] = magnitude_value;

        local_max_magnitude = std::max(magnitude_value, local_max_magnitude);
      }
    }

#pragma omp critical
    {
      max_magnitude = std::max(local_max_magnitude, max_magnitude);
    }
  }
}

void BadanovASelectEdgeSobelOMP::ComputeGradientAtPixel(const std::vector<uint8_t> &input, int row, int col,
                                                        float &gradient_x, float &gradient_y) const {
  gradient_x = 0.0F;
  gradient_y = 0.0F;

  for (int kernel_row = -1; kernel_row <= 1; ++kernel_row) {
    for (int kernel_col = -1; kernel_col <= 1; ++kernel_col) {
      const size_t pixel_index =
          (static_cast<size_t>(row + kernel_row) * static_cast<size_t>(width_)) + static_cast<size_t>(col + kernel_col);
      const uint8_t pixel = input[pixel_index];

      const int kx_idx = kernel_row + 1;
      const int ky_idx = kernel_col + 1;
      const int kernel_x_value = kKernelX.at(static_cast<size_t>(kx_idx)).at(static_cast<size_t>(ky_idx));
      const int kernel_y_value = kKernelY.at(static_cast<size_t>(kx_idx)).at(static_cast<size_t>(ky_idx));

      gradient_x += static_cast<float>(pixel) * static_cast<float>(kernel_x_value);
      gradient_y += static_cast<float>(pixel) * static_cast<float>(kernel_y_value);
    }
  }
}

void BadanovASelectEdgeSobelOMP::ApplyThreshold(const std::vector<float> &magnitude, float max_magnitude,
                                                std::vector<uint8_t> &output) const {
  if (max_magnitude > 0.0F) {
    const float scale = 255.0F / max_magnitude;
    const size_t size = magnitude.size();
    const int threshold = threshold_;

#pragma omp parallel for schedule(static) default(none) shared(magnitude, output, scale, size) firstprivate(threshold)
    for (size_t i = 0; i < size; ++i) {
      output[i] = (magnitude[i] * scale > static_cast<float>(threshold_)) ? 255 : 0;
    }
  } else {
    std::ranges::fill(output, 0);
  }
}

bool BadanovASelectEdgeSobelOMP::RunImpl() {
  const auto &input = GetInput();
  auto &output = GetOutput();

  if (height_ < 3 || width_ < 3) {
    output = input;
    return true;
  }

  std::vector<float> magnitude(input.size(), 0.0F);
  float max_magnitude = 0.0F;

  ApplySobelOperator(input, magnitude, max_magnitude);
  ApplyThreshold(magnitude, max_magnitude, output);

  return true;
}

bool BadanovASelectEdgeSobelOMP::PostProcessingImpl() {
  return true;
}

}  // namespace badanov_a_select_edge_sobel
