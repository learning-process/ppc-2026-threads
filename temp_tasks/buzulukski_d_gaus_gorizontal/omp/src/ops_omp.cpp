#include "buzulukski_d_gaus_gorizontal/omp/include/ops_omp.hpp"

#include <omp.h>

#include <algorithm>
#include <array>
#include <cstddef>
#include <cstdint>
#include <vector>

#include "buzulukski_d_gaus_gorizontal/common/include/common.hpp"

namespace buzulukski_d_gaus_gorizontal {

namespace {
constexpr int kChannels = 3;
constexpr int kKernelSize = 3;
constexpr int kKernelSum = 16;

using KernelRow = std::array<int, kKernelSize>;
constexpr std::array<KernelRow, kKernelSize> kKernel = {{{1, 2, 1}, {2, 4, 2}, {1, 2, 1}}};

uint8_t CalculatePixel(const uint8_t *in, int py, int px, int w, int h, int ch) {
  int sum = 0;
  for (int ky = -1; ky <= 1; ++ky) {
    for (int kx = -1; kx <= 1; ++kx) {
      int ny = std::max(0, std::min(h - 1, py + ky));
      int nx = std::max(0, std::min(w - 1, px + kx));

      size_t idx = (((static_cast<size_t>(ny) * static_cast<size_t>(w)) + static_cast<size_t>(nx)) * 3) +
                   static_cast<size_t>(ch);

      size_t row_idx = static_cast<size_t>(ky) + 1;
      size_t col_idx = static_cast<size_t>(kx) + 1;
      sum += static_cast<int>(in[idx]) * kKernel.at(row_idx).at(col_idx);
    }
  }
  return static_cast<uint8_t>(sum / kKernelSum);
}
}  // namespace

BuzulukskiDGausGorizontalOMP::BuzulukskiDGausGorizontalOMP(const InType &in) : BaseTask() {
  SetTypeOfTask(GetStaticTypeOfTask());
  GetInput() = in;
  GetOutput() = 0;
}

bool BuzulukskiDGausGorizontalOMP::ValidationImpl() {
  return GetInput() >= kKernelSize;
}

bool BuzulukskiDGausGorizontalOMP::PreProcessingImpl() {
  width_ = GetInput();
  height_ = GetInput();

  if (width_ < kKernelSize) {
    return false;
  }

  const auto total_size = (static_cast<std::size_t>(width_) * static_cast<std::size_t>(height_) * kChannels);
  input_image_.assign(total_size, static_cast<uint8_t>(100));
  output_image_.assign(total_size, 0);
  return true;
}

void BuzulukskiDGausGorizontalOMP::ApplyGaussianToPixel(int py, int px) {
  (void)py;
  (void)px;
}

bool BuzulukskiDGausGorizontalOMP::RunImpl() {
  const int h = height_;
  const int w = width_;
  const uint8_t *in_ptr = input_image_.data();
  uint8_t *out_ptr = output_image_.data();

#pragma omp parallel for default(none) shared(h, w, in_ptr, out_ptr, kKernel)
  for (int py = 0; py < h; ++py) {
    for (int px = 0; px < w; ++px) {
      for (int ch = 0; ch < 3; ++ch) {
        size_t out_idx = ((((static_cast<size_t>(py) * static_cast<size_t>(w)) + static_cast<size_t>(px)) * 3) +
                          static_cast<size_t>(ch));
        out_ptr[out_idx] = CalculatePixel(in_ptr, py, px, w, h, ch);
      }
    }
  }
  return true;
}

bool BuzulukskiDGausGorizontalOMP::PostProcessingImpl() {
  if (output_image_.empty()) {
    return false;
  }

  int64_t total_sum = 0;
  for (const auto &val : output_image_) {
    total_sum += static_cast<int64_t>(val);
  }
  GetOutput() = static_cast<int>(total_sum / static_cast<int64_t>(output_image_.size()));
  return true;
}

}  // namespace buzulukski_d_gaus_gorizontal
