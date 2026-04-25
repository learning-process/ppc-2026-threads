#include "vdovin_a_gauss_block/stl/include/ops_stl.hpp"

#include <algorithm>
#include <array>
#include <cstdint>
#include <thread>
#include <vector>

#include "util/include/util.hpp"
#include "vdovin_a_gauss_block/common/include/common.hpp"

namespace vdovin_a_gauss_block {

namespace {
constexpr int kChannels = 3;
constexpr int kKernelSize = 3;
constexpr int kKernelSum = 16;
constexpr std::array<std::array<int, kKernelSize>, kKernelSize> kKernel = {{{1, 2, 1}, {2, 4, 2}, {1, 2, 1}}};
}  // namespace

VdovinAGaussBlockSTL::VdovinAGaussBlockSTL(const InType &in) {
  SetTypeOfTask(GetStaticTypeOfTask());
  GetInput() = in;
  GetOutput() = 0;
}

bool VdovinAGaussBlockSTL::ValidationImpl() {
  return GetInput() >= 3;
}

bool VdovinAGaussBlockSTL::PreProcessingImpl() {
  width_ = GetInput();
  height_ = GetInput();
  if (width_ < 3 || height_ < 3) {
    input_image_.clear();
    output_image_.clear();
    return false;
  }
  int total = width_ * height_ * kChannels;
  input_image_.assign(total, 100);
  output_image_.assign(total, 0);
  return true;
}

void VdovinAGaussBlockSTL::ApplyGaussianToPixel(int py, int px) {
  for (int ch = 0; ch < kChannels; ch++) {
    int sum = 0;
    for (int ky = -1; ky <= 1; ky++) {
      for (int kx = -1; kx <= 1; kx++) {
        int ny = std::clamp(py + ky, 0, height_ - 1);
        int nx = std::clamp(px + kx, 0, width_ - 1);
        sum += input_image_[(((ny * width_) + nx) * kChannels) + ch] * kKernel.at(ky + 1).at(kx + 1);
      }
    }
    output_image_[(((py * width_) + px) * kChannels) + ch] = static_cast<uint8_t>(std::clamp(sum / kKernelSum, 0, 255));
  }
}

void VdovinAGaussBlockSTL::ProcessRows(int row_start, int row_end) {
  for (int py = row_start; py < row_end; py++) {
    for (int px = 0; px < width_; px++) {
      ApplyGaussianToPixel(py, px);
    }
  }
}

bool VdovinAGaussBlockSTL::RunImpl() {
  if (input_image_.empty() || output_image_.empty()) {
    return false;
  }

  int num_threads = ppc::util::GetNumThreads();
  int rows_per_thread = height_ / num_threads;
  int remainder = height_ % num_threads;

  std::vector<std::thread> threads;
  threads.reserve(num_threads);

  int offset = 0;
  for (int ti = 0; ti < num_threads; ti++) {
    int count = rows_per_thread + (ti < remainder ? 1 : 0);
    threads.emplace_back(&VdovinAGaussBlockSTL::ProcessRows, this, offset, offset + count);
    offset += count;
  }

  for (auto &th : threads) {
    th.join();
  }

  return true;
}

bool VdovinAGaussBlockSTL::PostProcessingImpl() {
  if (output_image_.empty()) {
    return false;
  }
  auto total = static_cast<int64_t>(output_image_.size());
  if (total == 0) {
    return false;
  }
  int64_t sum = 0;
  for (int64_t idx = 0; idx < total; idx++) {
    sum += output_image_[idx];
  }
  GetOutput() = static_cast<int>(sum / total);
  return true;
}

}  // namespace vdovin_a_gauss_block
