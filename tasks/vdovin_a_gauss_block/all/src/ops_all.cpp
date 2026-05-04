#include "vdovin_a_gauss_block/all/include/ops_all.hpp"

#include <mpi.h>

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
constexpr std::array<std::array<int, kKernelSize>, kKernelSize> kKernel = {{{1, 2, 1}, {2, 4, 2}, {1, 2, 1}}};
}  // namespace

VdovinAGaussBlockALL::VdovinAGaussBlockALL(const InType &in) {
  SetTypeOfTask(GetStaticTypeOfTask());
  GetInput() = in;
  GetOutput() = 0;
}

bool VdovinAGaussBlockALL::ValidationImpl() {
  return GetInput() >= 3;
}

bool VdovinAGaussBlockALL::PreProcessingImpl() {
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

uint8_t VdovinAGaussBlockALL::ComputePixelChannel(const std::vector<uint8_t> &img, int width, int height, int py,
                                                  int px, int ch) {
  int sum = 0;
  for (int ky = -1; ky <= 1; ky++) {
    for (int kx = -1; kx <= 1; kx++) {
      int ny = std::clamp(py + ky, 0, height - 1);
      int nx = std::clamp(px + kx, 0, width - 1);
      sum += img[(((ny * width) + nx) * kChannels) + ch] * kKernel.at(ky + 1).at(kx + 1);
    }
  }
  return static_cast<uint8_t>(std::clamp(sum / kKernelSum, 0, 255));
}

bool VdovinAGaussBlockALL::RunImpl() {
  if (input_image_.empty() || output_image_.empty()) {
    return false;
  }

  int rank = 0;
  int size = 1;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &size);

  MPI_Bcast(&width_, 1, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Bcast(&height_, 1, MPI_INT, 0, MPI_COMM_WORLD);

  int total_pixels = width_ * height_;
  MPI_Bcast(input_image_.data(), total_pixels * kChannels, MPI_UNSIGNED_CHAR, 0, MPI_COMM_WORLD);

  int rows_per_proc = height_ / size;
  int extra = height_ % size;

  std::vector<int> send_counts(size);
  std::vector<int> displs(size);
  int offset = 0;
  for (int pi = 0; pi < size; pi++) {
    int rows = rows_per_proc + (pi < extra ? 1 : 0);
    send_counts[pi] = rows * width_ * kChannels;
    displs[pi] = offset;
    offset += send_counts[pi];
  }

  int local_rows = rows_per_proc + (rank < extra ? 1 : 0);
  int local_size = local_rows * width_ * kChannels;
  std::vector<uint8_t> local_output(local_size, 0);

  int row_start = 0;
  for (int pi = 0; pi < rank; pi++) {
    row_start += rows_per_proc + (pi < extra ? 1 : 0);
  }

#pragma omp parallel for schedule(static) default(none) shared(local_output, local_rows, row_start)
  for (int ry = 0; ry < local_rows; ry++) {
    int py = row_start + ry;
    for (int px = 0; px < width_; px++) {
      for (int ch = 0; ch < kChannels; ch++) {
        local_output[(((ry * width_) + px) * kChannels) + ch] =
            ComputePixelChannel(input_image_, width_, height_, py, px, ch);
      }
    }
  }

  MPI_Allgatherv(local_output.data(), local_size, MPI_UNSIGNED_CHAR, output_image_.data(), send_counts.data(),
                 displs.data(), MPI_UNSIGNED_CHAR, MPI_COMM_WORLD);

  return true;
}

bool VdovinAGaussBlockALL::PostProcessingImpl() {
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
