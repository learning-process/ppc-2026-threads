#include "badanov_a_select_edge_sobel/tbb/include/ops_tbb.hpp"

#include <tbb/blocked_range.h>
#include <tbb/parallel_for.h>
#include <tbb/parallel_reduce.h>

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <vector>

#include "badanov_a_select_edge_sobel/common/include/common.hpp"

namespace badanov_a_select_edge_sobel {

BadanovASelectEdgeSobelTBB::BadanovASelectEdgeSobelTBB(const InType &in) {
  SetTypeOfTask(GetStaticTypeOfTask());
  GetInput() = in;
  GetOutput() = std::vector<uint8_t>();
}

bool BadanovASelectEdgeSobelTBB::ValidationImpl() {
  const auto &input = GetInput();
  return !input.empty();
}

bool BadanovASelectEdgeSobelTBB::PreProcessingImpl() {
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

void BadanovASelectEdgeSobelTBB::ApplySobelOperator(const std::vector<uint8_t> &input, std::vector<float> &magnitude,
                                                    float &max_magnitude) {
  const int height = height_;
  const int width = width_;

  using Range = tbb::blocked_range<int>;

  max_magnitude = tbb::parallel_reduce(Range(1, height - 1), 0.0F, [&](const Range &r, float init) -> float {
    float local_max = init;
    for (int row = r.begin(); row < r.end(); ++row) {
      for (int col = 1; col < width - 1; ++col) {
        float gradient_x = 0.0F;
        float gradient_y = 0.0F;

        ComputeGradientAtPixel(input, row, col, gradient_x, gradient_y);

        const float magnitude_value = std::sqrt((gradient_x * gradient_x) + (gradient_y * gradient_y));
        const size_t idx = (static_cast<size_t>(row) * static_cast<size_t>(width)) + static_cast<size_t>(col);
        magnitude[idx] = magnitude_value;

        if (magnitude_value > local_max) {
          local_max = magnitude_value;
        }
      }
    }
    return local_max;
  }, [](float a, float b) -> float { return std::max(a, b); });
}

void BadanovASelectEdgeSobelTBB::ComputeGradientAtPixel(const std::vector<uint8_t> &input, int row, int col,
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

void BadanovASelectEdgeSobelTBB::ApplyThreshold(const std::vector<float> &magnitude, float max_magnitude,
                                                std::vector<uint8_t> &output) const {
  if (max_magnitude > 0.0F) {
    const float scale = 255.0F / max_magnitude;
    const size_t size = magnitude.size();
    const int threshold = threshold_;

    tbb::parallel_for(tbb::blocked_range<size_t>(0, size), [&](const tbb::blocked_range<size_t> &r) {
      for (size_t i = r.begin(); i < r.end(); ++i) {
        output[i] = (magnitude[i] * scale > static_cast<float>(threshold)) ? 255 : 0;
      }
    });
  } else {
    std::ranges::fill(output, 0);
  }
}

bool BadanovASelectEdgeSobelTBB::RunImpl() {
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

bool BadanovASelectEdgeSobelTBB::PostProcessingImpl() {
  return true;
}

}  // namespace badanov_a_select_edge_sobel
