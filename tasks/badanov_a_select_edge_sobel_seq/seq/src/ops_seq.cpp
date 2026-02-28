#include "badanov_a_select_edge_sobel_seq/seq/include/ops_seq.hpp"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <vector>

#include "badanov_a_select_edge_sobel_seq/common/include/common.hpp"

namespace badanov_a_select_edge_sobel_seq {

BadanovASelectEdgeSobelSEQ::BadanovASelectEdgeSobelSEQ(const InType &in) {
  SetTypeOfTask(GetStaticTypeOfTask());
  GetInput() = in;
  GetOutput() = std::vector<uint8_t>();
}

bool BadanovASelectEdgeSobelSEQ::ValidationImpl() {
  const auto &input = GetInput();
  return !input.empty();
}

bool BadanovASelectEdgeSobelSEQ::PreProcessingImpl() {
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

void BadanovASelectEdgeSobelSEQ::ApplySobelOperator(const std::vector<uint8_t> &input, std::vector<float> &magnitude,
                                                    float &max_magnitude) {
  max_magnitude = 0.0F;

  for (int row = 1; row < height_ - 1; ++row) {
    for (int col = 1; col < width_ - 1; ++col) {
      float gradient_x = 0.0F;
      float gradient_y = 0.0F;

      ComputeGradientAtPixel(input, row, col, gradient_x, gradient_y);

      const float magnitude_value = std::sqrt((gradient_x * gradient_x) + (gradient_y * gradient_y));
      const size_t idx = (static_cast<size_t>(row) * static_cast<size_t>(width_)) + static_cast<size_t>(col);
      magnitude[idx] = magnitude_value;

      max_magnitude = std::max(magnitude_value, max_magnitude);
    }
  }
}

void BadanovASelectEdgeSobelSEQ::ComputeGradientAtPixel(const std::vector<uint8_t> &input, int row, int col,
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

void BadanovASelectEdgeSobelSEQ::ApplyThreshold(const std::vector<float> &magnitude, float max_magnitude,
                                                std::vector<uint8_t> &output) const {
  if (max_magnitude > 0.0F) {
    const float scale = 255.0F / max_magnitude;
    for (size_t i = 0; i < magnitude.size(); ++i) {
      output[i] = (magnitude[i] * scale > static_cast<float>(threshold_)) ? 255 : 0;
    }
  } else {
    std::ranges::fill(output, 0);
  }
}

bool BadanovASelectEdgeSobelSEQ::RunImpl() {
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

bool BadanovASelectEdgeSobelSEQ::PostProcessingImpl() {
  return true;
}

}  // namespace badanov_a_select_edge_sobel_seq
