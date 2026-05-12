#include "buzulukski_d_gaus_gorizontal/tbb/include/ops_tbb.hpp"

#include <tbb/blocked_range.h>
#include <tbb/parallel_for.h>

#include <algorithm>
#include <array>
#include <cstddef>
#include <cstdint>

#include "buzulukski_d_gaus_gorizontal/common/include/common.hpp"

namespace buzulukski_d_gaus_gorizontal {

namespace {
constexpr int kChannels = 3;
constexpr int kKernelSize = 3;
constexpr int kKernelSum = 16;

using KernelRow = std::array<int, kKernelSize>;
constexpr std::array<KernelRow, kKernelSize> kKernel = {{{1, 2, 1}, {2, 4, 2}, {1, 2, 1}}};

uint8_t CalculatePixelTBB(const uint8_t *in, int py, int px, int w, int h, int ch) {
  int sum = 0;
  for (int ky = -1; ky <= 1; ++ky) {
    for (int kx = -1; kx <= 1; ++kx) {
      int ny = std::clamp(py + ky, 0, h - 1);
      int nx = std::clamp(px + kx, 0, w - 1);

      size_t idx = (((static_cast<size_t>(ny) * static_cast<size_t>(w)) + static_cast<size_t>(nx)) * kChannels) +
                   static_cast<size_t>(ch);

      size_t row_idx = static_cast<size_t>(ky) + 1;
      size_t col_idx = static_cast<size_t>(kx) + 1;
      sum += static_cast<int>(in[idx]) * kKernel.at(row_idx).at(col_idx);
    }
  }
  return static_cast<uint8_t>(sum / kKernelSum);
}
}  // namespace

BuzulukskiDGausGorizontalTBB::BuzulukskiDGausGorizontalTBB(const InType &in) : BaseTask() {
  SetTypeOfTask(GetStaticTypeOfTask());
  GetInput() = in;
  GetOutput() = 0;
}

bool BuzulukskiDGausGorizontalTBB::ValidationImpl() {
  return GetInput() >= kKernelSize;
}

bool BuzulukskiDGausGorizontalTBB::PreProcessingImpl() {
  width_ = GetInput();
  height_ = GetInput();
  const auto total_size = static_cast<std::size_t>(width_) * static_cast<std::size_t>(height_) * kChannels;
  input_image_.assign(total_size, 100);
  output_image_.assign(total_size, 0);
  return true;
}

bool BuzulukskiDGausGorizontalTBB::RunImpl() {
  const int h = height_;
  const int w = width_;
  const uint8_t *in_ptr = input_image_.data();
  uint8_t *out_ptr = output_image_.data();

  tbb::parallel_for(tbb::blocked_range<int>(0, h), [&](const tbb::blocked_range<int> &r) {
    for (int py = r.begin(); py < r.end(); ++py) {
      for (int px = 0; px < w; ++px) {
        for (int ch = 0; ch < kChannels; ++ch) {
          size_t out_idx =
              (((static_cast<size_t>(py) * static_cast<size_t>(w)) + static_cast<size_t>(px)) * kChannels) +
              static_cast<size_t>(ch);
          out_ptr[out_idx] = CalculatePixelTBB(in_ptr, py, px, w, h, ch);
        }
      }
    }
  });
  return true;
}

bool BuzulukskiDGausGorizontalTBB::PostProcessingImpl() {
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
