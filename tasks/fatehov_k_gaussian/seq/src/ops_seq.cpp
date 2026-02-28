#include "fatehov_k_gaussian/seq/include/ops_seq.hpp"

#include <algorithm>
#include <cmath>
#include <numeric>

#include "util/include/util.hpp"

namespace fatehov_k_gaussian {

FatehovKGaussianSEQ::FatehovKGaussianSEQ(const InType &in) {
  SetTypeOfTask(GetStaticTypeOfTask());
  GetInput() = in;
}

bool FatehovKGaussianSEQ::ValidationImpl() {
  const auto &input = GetInput();
  return input.image.width > 0 && input.image.height > 0 && input.image.channels > 0 && !input.image.data.empty() &&
         input.sigma > 0.0f;
}

bool FatehovKGaussianSEQ::PreProcessingImpl() {
  const auto &input = GetInput();
  float sigma = input.sigma;

  // Определение размера ядра по правилу 3-х сигм (2 * ceil(3*sigma) + 1)
  kernel_size_ = 2 * static_cast<int>(std::ceil(3.0f * sigma)) + 1;
  kernel_.resize(static_cast<size_t>(kernel_size_) * kernel_size_);

  int half = kernel_size_ / 2;
  float sum = 0.0f;
  float two_sigma_sq = 2.0f * sigma * sigma;

  for (int i = -half; i <= half; i++) {
    for (int j = -half; j <= half; j++) {
      float val = std::exp(-(static_cast<float>(i * i + j * j)) / two_sigma_sq);
      kernel_[(i + half) * kernel_size_ + (j + half)] = val;
      sum += val;
    }
  }

  // Нормализация, чтобы сумма коэффициентов была равна 1
  for (float &val : kernel_) {
    val /= sum;
  }

  auto &output = GetOutput();
  output = Image(input.image.width, input.image.height, input.image.channels);

  return true;
}

bool FatehovKGaussianSEQ::RunImpl() {
  const auto &input = GetInput();
  auto &output = GetOutput();
  const int w = static_cast<int>(input.image.width);
  const int h = static_cast<int>(input.image.height);
  const int ch = static_cast<int>(input.image.channels);
  const int half = kernel_size_ / 2;

  for (int y = 0; y < h; y++) {
    for (int x = 0; x < w; x++) {
      for (int c = 0; c < ch; c++) {
        float res = 0.0f;
        for (int ky = -half; ky <= half; ky++) {
          for (int kx = -half; kx <= half; kx++) {
            // Метод Clamp для обработки границ
            int ny = std::clamp(y + ky, 0, h - 1);
            int nx = std::clamp(x + kx, 0, w - 1);

            float weight = kernel_[(ky + half) * kernel_size_ + (kx + half)];
            res += static_cast<float>(input.image.data[(ny * w + nx) * ch + c]) * weight;
          }
        }
        output.data[(y * w + x) * ch + c] = static_cast<uint8_t>(std::clamp(res, 0.0f, 255.0f));
      }
    }
  }
  return true;
}

bool FatehovKGaussianSEQ::PostProcessingImpl() {
  return true;
}

}  // namespace fatehov_k_gaussian
