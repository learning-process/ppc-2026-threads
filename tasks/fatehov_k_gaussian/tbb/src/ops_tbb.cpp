#include "fatehov_k_gaussian/tbb/include/ops_tbb.hpp"

#include <oneapi/tbb/blocked_range.h>
#include <oneapi/tbb/parallel_for.h>

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <cstdint>

namespace fatehov_k_gaussian {

FatehovKGaussianTBB::FatehovKGaussianTBB(const InType &in) {
  SetTypeOfTask(GetStaticTypeOfTask());
  GetInput() = in;
  GetOutput() = Image{};
}

bool FatehovKGaussianTBB::ValidationImpl() {
  const auto &input = GetInput();

  if (input.sigma <= 0.0F) {
    return false;
  }

  if (input.image.width == 0 || input.image.height == 0 || input.image.channels == 0) {
    return false;
  }

  const std::size_t expected_size =
      static_cast<std::size_t>(input.image.width) * input.image.height * input.image.channels;

  return input.image.data.size() == expected_size;
}

bool FatehovKGaussianTBB::PreProcessingImpl() {
  const auto &input_image = GetInput().image;

  GetOutput() = Image(input_image.width, input_image.height, input_image.channels);

  return GetOutput().data.size() == input_image.data.size();
}

bool FatehovKGaussianTBB::RunImpl() {
  const auto &input = GetInput();
  const auto &src = input.image;
  auto &dst = GetOutput();

  const uint32_t width = src.width;
  const uint32_t height = src.height;
  const uint32_t channels = src.channels;

  const float sigma = input.sigma;

  std::array<float, 9> kernel{};
  float kernel_sum = 0.0F;

  int idx = 0;
  for (int y = -1; y <= 1; ++y) {
    for (int x = -1; x <= 1; ++x) {
      const float value = std::exp(-(static_cast<float>(x * x + y * y)) / (2.0F * sigma * sigma));
      kernel[idx++] = value;
      kernel_sum += value;
    }
  }

  for (auto &value : kernel) {
    value /= kernel_sum;
  }

  oneapi::tbb::parallel_for(oneapi::tbb::blocked_range<std::size_t>(0, static_cast<std::size_t>(height)),
                            [&](const oneapi::tbb::blocked_range<std::size_t> &range) {
    for (std::size_t y = range.begin(); y < range.end(); ++y) {
      for (std::size_t x = 0; x < width; ++x) {
        for (std::size_t c = 0; c < channels; ++c) {
          float pixel_value = 0.0F;

          int kernel_index = 0;
          for (int ky = -1; ky <= 1; ++ky) {
            const int current_y = std::clamp(static_cast<int>(y) + ky, 0, static_cast<int>(height) - 1);

            for (int kx = -1; kx <= 1; ++kx) {
              const int current_x = std::clamp(static_cast<int>(x) + kx, 0, static_cast<int>(width) - 1);

              const std::size_t src_index = (static_cast<std::size_t>(current_y) * width + current_x) * channels + c;

              pixel_value += static_cast<float>(src.data[src_index]) * kernel[kernel_index];
              ++kernel_index;
            }
          }

          const std::size_t dst_index = (y * width + x) * channels + c;

          const int rounded_value = static_cast<int>(std::round(pixel_value));
          dst.data[dst_index] = static_cast<uint8_t>(std::clamp(rounded_value, 0, 255));
        }
      }
    }
  });

  return true;
}

bool FatehovKGaussianTBB::PostProcessingImpl() {
  return !GetOutput().data.empty();
}

}  // namespace fatehov_k_gaussian
