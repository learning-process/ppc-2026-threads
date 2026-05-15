#include "fatehov_k_gaussian/tbb/include/ops_tbb.hpp"

#include <oneapi/tbb/blocked_range.h>
#include <oneapi/tbb/parallel_for.h>

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <vector>

#include "fatehov_k_gaussian/common/include/common.hpp"

namespace fatehov_k_gaussian {

FatehovKGaussianTBB::FatehovKGaussianTBB(const InType &in) {
  SetTypeOfTask(GetStaticTypeOfTask());
  GetInput() = in;
}

bool FatehovKGaussianTBB::ValidationImpl() {
  const auto &input = GetInput();
  return input.image.width > 0 && input.image.height > 0 && input.image.channels > 0 && !input.image.data.empty() &&
         input.sigma > 0.0F;
}

bool FatehovKGaussianTBB::PreProcessingImpl() {
  const auto &input = GetInput();
  const float sigma = input.sigma;

  kernel_size_ = (2 * static_cast<int>(std::ceil(3.0F * sigma))) + 1;
  kernel_.resize(static_cast<std::size_t>(kernel_size_) * kernel_size_);

  const int half = kernel_size_ / 2;
  const float two_sigma_sq = 2.0F * sigma * sigma;
  float sum = 0.0F;

  for (int i = -half; i <= half; ++i) {
    for (int j = -half; j <= half; ++j) {
      const float val = std::exp(-(static_cast<float>((i * i) + (j * j))) / two_sigma_sq);
      kernel_[((i + half) * kernel_size_) + (j + half)] = val;
      sum += val;
    }
  }

  for (float &val : kernel_) {
    val /= sum;
  }

  GetOutput() = Image(input.image.width, input.image.height, input.image.channels);

  return true;
}

bool FatehovKGaussianTBB::RunImpl() {
  const auto &input = GetInput();
  auto &output = GetOutput();
  const int w = static_cast<int>(input.image.width);
  const int h = static_cast<int>(input.image.height);
  const int ch = static_cast<int>(input.image.channels);
  const int half = kernel_size_ / 2;
  const int kernel_size = kernel_size_;
  const auto &kernel = kernel_;

  oneapi::tbb::parallel_for(oneapi::tbb::blocked_range<int>(0, h), [&](const oneapi::tbb::blocked_range<int> &range) {
    for (int y_coord = range.begin(); y_coord < range.end(); ++y_coord) {
      for (int x_coord = 0; x_coord < w; ++x_coord) {
        for (int c_coord = 0; c_coord < ch; ++c_coord) {
          float res = 0.0F;
          for (int ky = -half; ky <= half; ++ky) {
            for (int kx = -half; kx <= half; ++kx) {
              const int ny = std::clamp(y_coord + ky, 0, h - 1);
              const int nx = std::clamp(x_coord + kx, 0, w - 1);
              const float weight = kernel[((ky + half) * kernel_size) + (kx + half)];
              res += static_cast<float>(input.image.data[((ny * w + nx) * ch) + c_coord]) * weight;
            }
          }
          output.data[((y_coord * w + x_coord) * ch) + c_coord] = static_cast<uint8_t>(std::clamp(res, 0.0F, 255.0F));
        }
      }
    }
  });

  return true;
}

bool FatehovKGaussianTBB::PostProcessingImpl() {
  return true;
}

}  // namespace fatehov_k_gaussian
