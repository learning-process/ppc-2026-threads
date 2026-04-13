#include "vdovin_a_gauss_block/tbb/include/ops_tbb.hpp"

#include <algorithm>
#include <array>
#include <cstdint>
#include <vector>

#include "oneapi/tbb/blocked_range2d.h"
#include "oneapi/tbb/parallel_for.h"
#include "vdovin_a_gauss_block/common/include/common.hpp"

namespace vdovin_a_gauss_block {

namespace {
constexpr int kChannels = 3;
constexpr int kKernelSize = 3;
constexpr int kKernelSum = 16;
constexpr std::array<std::array<int, kKernelSize>, kKernelSize> kKernel = {{{1, 2, 1}, {2, 4, 2}, {1, 2, 1}}};
}  // namespace

VdovinAGaussBlockTBB::VdovinAGaussBlockTBB(const InType &in) {
  SetTypeOfTask(GetStaticTypeOfTask());
  GetInput() = in;
  GetOutput() = 0;
}

bool VdovinAGaussBlockTBB::ValidationImpl() {
  return GetInput() >= 3;
}

bool VdovinAGaussBlockTBB::PreProcessingImpl() {
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

void VdovinAGaussBlockTBB::ApplyGaussianToPixel(int py, int px) {
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

bool VdovinAGaussBlockTBB::RunImpl() {
  if (input_image_.empty() || output_image_.empty()) {
    return false;
  }

  int block_height = std::max(1, height_ / 4);
  int block_width = std::max(1, width_ / 4);

  tbb::parallel_for(tbb::blocked_range2d<int>(0, height_, block_height, 0, width_, block_width),
                    [this](const tbb::blocked_range2d<int> &range) {
    for (int py = range.rows().begin(); py < range.rows().end(); py++) {
      for (int px = range.cols().begin(); px < range.cols().end(); px++) {
        ApplyGaussianToPixel(py, px);
      }
    }
  });

  return true;
}

bool VdovinAGaussBlockTBB::PostProcessingImpl() {
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
