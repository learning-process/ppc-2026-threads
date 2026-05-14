#include "buzulukski_d_gaus_gorizontal/stl/include/ops_stl.hpp"

#include <algorithm>
#include <array>
#include <cstddef>
#include <cstdint>
#include <thread>
#include <vector>

#include "buzulukski_d_gaus_gorizontal/common/include/common.hpp"

namespace buzulukski_d_gaus_gorizontal {

namespace {
constexpr int kChannels = 3;
constexpr int kKernelSize = 3;
constexpr int kKernelSum = 16;

using KernelRow = std::array<int, kKernelSize>;
constexpr std::array<KernelRow, kKernelSize> kKernel = {{{1, 2, 1}, {2, 4, 2}, {1, 2, 1}}};

uint8_t CalculatePixelSTL(const uint8_t *in, int py, int px, int w, int h, int ch) {
  int sum = 0;
  for (int ky = -1; ky <= 1; ++ky) {
    for (int kx = -1; kx <= 1; ++kx) {
      int ny = std::clamp(py + ky, 0, h - 1);
      int nx = std::clamp(px + kx, 0, w - 1);

      size_t idx = ((static_cast<size_t>(ny) * w + nx) * kChannels) + ch;
      sum += static_cast<int>(in[idx]) * kKernel.at(ky + 1).at(kx + 1);
    }
  }
  return static_cast<uint8_t>(sum / kKernelSum);
}

void ProcessRows(int start_row, int end_row, int w, int h, const uint8_t *in, uint8_t *out) {
  for (int py = start_row; py < end_row; ++py) {
    for (int px = 0; px < w; ++px) {
      for (int ch = 0; ch < kChannels; ++ch) {
        size_t out_idx = ((static_cast<size_t>(py) * w + px) * kChannels) + ch;
        out[out_idx] = CalculatePixelSTL(in, py, px, w, h, ch);
      }
    }
  }
}
}  // namespace

BuzulukskiDGausGorizontalSTL::BuzulukskiDGausGorizontalSTL(const InType &in) : BaseTask() {
  SetTypeOfTask(GetStaticTypeOfTask());
  GetInput() = in;
  GetOutput() = 0;
}

bool BuzulukskiDGausGorizontalSTL::ValidationImpl() {
  return GetInput() >= kKernelSize;
}

bool BuzulukskiDGausGorizontalSTL::PreProcessingImpl() {
  width_ = GetInput();
  height_ = GetInput();
  const auto total_size = static_cast<std::size_t>(width_) * height_ * kChannels;
  input_image_.assign(total_size, 100);
  output_image_.assign(total_size, 0);
  return true;
}

bool BuzulukskiDGausGorizontalSTL::RunImpl() {
  const int h = height_;
  const int w = width_;
  const uint8_t *in_ptr = input_image_.data();
  uint8_t *out_ptr = output_image_.data();

  unsigned int raw_threads = std::thread::hardware_concurrency();
  int n_threads = (raw_threads == 0) ? 2 : static_cast<int>(raw_threads);

  std::vector<std::thread> threads;
  threads.reserve(n_threads);

  int rows_per_thread = h / n_threads;

  for (int i = 0; i < n_threads; ++i) {
    int start = i * rows_per_thread;
    int end = (i == n_threads - 1) ? h : (i + 1) * rows_per_thread;

    threads.emplace_back(ProcessRows, start, end, w, h, in_ptr, out_ptr);
  }

  for (auto &t : threads) {
    if (t.joinable()) {
      t.join();
    }
  }

  return true;
}

bool BuzulukskiDGausGorizontalSTL::PostProcessingImpl() {
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
