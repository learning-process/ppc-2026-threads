#include "vdovin_a_gauss_block/omp/include/ops_omp.hpp"

#include <algorithm>
#include <array>
#include <cstdint>
#include <vector>

#include "vdovin_a_gauss_block/common/include/common.hpp"

namespace vdovin_a_gauss_block {

namespace {
constexpr int kChannels = 3;
constexpr int kKernelSize = 3;
constexpr int kKernelSum = 16;
constexpr int kBlockSize = 32;
constexpr std::array<std::array<int, kKernelSize>, kKernelSize> kKernel = {{{1, 2, 1}, {2, 4, 2}, {1, 2, 1}}};
}  // namespace

VdovinAGaussBlockOMP::VdovinAGaussBlockOMP(const InType &in) {
  SetTypeOfTask(GetStaticTypeOfTask());
  GetInput() = in;
  GetOutput() = 0;
}

bool VdovinAGaussBlockOMP::ValidationImpl() {
  return GetInput() >= 3;
}

bool VdovinAGaussBlockOMP::PreProcessingImpl() {
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

void VdovinAGaussBlockOMP::ApplyGaussianToPixel(int py, int px) {
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

bool VdovinAGaussBlockOMP::RunImpl() {
  if (input_image_.empty() || output_image_.empty()) {
    return false;
  }

  int num_blocks_y = (height_ + kBlockSize - 1) / kBlockSize;
  int num_blocks_x = (width_ + kBlockSize - 1) / kBlockSize;
  int total_blocks = num_blocks_y * num_blocks_x;

#pragma omp parallel for schedule(static) default(none) shared(total_blocks, num_blocks_x)
  for (int bi = 0; bi < total_blocks; bi++) {
    int by = (bi / num_blocks_x) * kBlockSize;
    int bx = (bi % num_blocks_x) * kBlockSize;
    int y_end = std::min(by + kBlockSize, height_);
    int x_end = std::min(bx + kBlockSize, width_);
    for (int py = by; py < y_end; py++) {
      for (int px = bx; px < x_end; px++) {
        ApplyGaussianToPixel(py, px);
      }
    }
  }

  return true;
}

bool VdovinAGaussBlockOMP::PostProcessingImpl() {
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
