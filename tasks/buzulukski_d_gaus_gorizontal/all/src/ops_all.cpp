#include "buzulukski_d_gaus_gorizontal/all/include/ops_all.hpp"

#include <mpi.h>

#include <algorithm>
#include <vector>

#include "buzulukski_d_gaus_gorizontal/common/include/common.hpp"

namespace buzulukski_d_gaus_gorizontal {

namespace {
constexpr int kChannels = 3;
constexpr int kKernelSum = 16;
const int kKernel[3][3] = {{1, 2, 1}, {2, 4, 2}, {1, 2, 1}};

uint8_t CalculatePixelALL(const uint8_t *in, int py, int px, int w, int h, int ch) {
  int sum = 0;
  for (int ky = -1; ky <= 1; ++ky) {
    for (int kx = -1; kx <= 1; ++kx) {
      int ny = std::clamp(py + ky, 0, h - 1);
      int nx = std::clamp(px + kx, 0, w - 1);
      size_t idx = (static_cast<size_t>(ny) * w + nx) * kChannels + ch;
      sum += static_cast<int>(in[idx]) * kKernel[ky + 1][kx + 1];
    }
  }
  return static_cast<uint8_t>(sum / kKernelSum);
}
}  // namespace

BuzulukskiDGausGorizontalALL::BuzulukskiDGausGorizontalALL(const InType &in) : BaseTask() {
  SetTypeOfTask(GetStaticTypeOfTask());
  GetInput() = in;
  GetOutput() = 0;
}

bool BuzulukskiDGausGorizontalALL::ValidationImpl() {
  return GetInput() >= 3;
}

bool BuzulukskiDGausGorizontalALL::PreProcessingImpl() {
  width_ = GetInput();
  height_ = GetInput();
  size_t total_size = static_cast<size_t>(width_) * height_ * kChannels;
  input_image_.assign(total_size, 100);
  output_image_.assign(total_size, 0);
  return true;
}

bool BuzulukskiDGausGorizontalALL::RunImpl() {
  int rank;
  int size;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &size);

  int h = height_;
  int w = width_;

  int rows_per_proc = h / size;
  int start_row = rank * rows_per_proc;
  int end_row = (rank == size - 1) ? h : (rank + 1) * rows_per_proc;

  std::vector<uint8_t> local_output((static_cast<size_t>(end_row) - start_row) * w * kChannels);

#pragma omp parallel for collapse(2)
  for (int py = start_row; py < end_row; ++py) {
    for (int px = 0; px < w; ++px) {
      for (int ch = 0; ch < kChannels; ++ch) {
        int local_py = py - start_row;
        size_t local_idx = (static_cast<size_t>(local_py) * w + px) * kChannels + ch;
        local_output[local_idx] = CalculatePixelALL(input_image_.data(), py, px, w, h, ch);
      }
    }
  }

  std::vector<int> recv_counts(size);
  std::vector<int> displs(size);

  if (rank == 0) {
    for (int i = 0; i < size; ++i) {
      int s = i * rows_per_proc;
      int e = (i == size - 1) ? h : (i + 1) * rows_per_proc;
      recv_counts[i] = (e - s) * w * kChannels;
      displs[i] = s * w * kChannels;
    }
  }

  MPI_Gatherv(local_output.data(), static_cast<int>(local_output.size()), MPI_UNSIGNED_CHAR, output_image_.data(),
              recv_counts.data(), displs.data(), MPI_UNSIGNED_CHAR, 0, MPI_COMM_WORLD);

  return true;
}

bool BuzulukskiDGausGorizontalALL::PostProcessingImpl() {
  int rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  if (rank == 0) {
    int64_t total_sum = 0;
    for (const auto &val : output_image_) {
      total_sum += val;
    }
    GetOutput() = static_cast<int>(total_sum / static_cast<int64_t>(output_image_.size()));
  }
  return true;
}

}  // namespace buzulukski_d_gaus_gorizontal
