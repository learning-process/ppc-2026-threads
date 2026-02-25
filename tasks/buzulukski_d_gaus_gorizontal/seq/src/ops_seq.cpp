#include "buzulukski_d_gaus_gorizontal/seq/include/ops_seq.hpp"

#include <algorithm>
#include <array>

namespace buzulukski_d_gaus_gorizontal {

namespace {
constexpr int kChannels = 3;
constexpr int kKernelSize = 3;
constexpr int kKernelSum = 16;
constexpr std::array<std::array<int, kKernelSize>, kKernelSize> kKernel = {{{1, 2, 1}, {2, 4, 2}, {1, 2, 1}}};
}  // namespace

// ИСПРАВЛЕНИЕ: Вызываем конструктор по умолчанию BaseTask()
BuzulukskiDGausGorizontalSEQ::BuzulukskiDGausGorizontalSEQ(const InType &in) : BaseTask() {
  SetTypeOfTask(GetStaticTypeOfTask());
  GetInput() = in;
  GetOutput() = 0;
}

bool BuzulukskiDGausGorizontalSEQ::ValidationImpl() {
  return GetInput() >= 3;
}

bool BuzulukskiDGausGorizontalSEQ::PreProcessingImpl() {
  width_ = GetInput();
  height_ = GetInput();
  if (width_ < 3) {
    return false;
  }

  int total = width_ * height_ * kChannels;
  input_image_.assign(total, 100);
  output_image_.assign(total, 0);
  return true;
}

void BuzulukskiDGausGorizontalSEQ::ApplyGaussianToPixel(int py, int px) {
  for (int ch = 0; ch < kChannels; ch++) {
    int sum = 0;
    for (int ky = -1; ky <= 1; ky++) {
      for (int kx = -1; kx <= 1; kx++) {
        int ny = std::clamp(py + ky, 0, height_ - 1);
        int nx = std::clamp(px + kx, 0, width_ - 1);
        sum += input_image_[(((ny * width_) + nx) * kChannels) + ch] * kKernel[ky + 1][kx + 1];
      }
    }
    output_image_[(((py * width_) + px) * kChannels) + ch] = static_cast<uint8_t>(sum / kKernelSum);
  }
}

bool BuzulukskiDGausGorizontalSEQ::RunImpl() {
  for (int py = 0; py < height_; py++) {
    for (int px = 0; px < width_; px++) {
      ApplyGaussianToPixel(py, px);
    }
  }
  return true;
}

bool BuzulukskiDGausGorizontalSEQ::PostProcessingImpl() {
  if (output_image_.empty()) {
    return false;
  }
  int64_t sum = 0;
  for (uint8_t val : output_image_) {
    sum += val;
  }
  GetOutput() = static_cast<int>(sum / output_image_.size());
  return true;
}

}  // namespace buzulukski_d_gaus_gorizontal
