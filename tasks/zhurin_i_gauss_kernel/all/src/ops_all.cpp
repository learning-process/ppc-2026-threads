#include "zhurin_i_gauss_kernel/all/include/ops_all.hpp"

#include <mpi.h>
#include <omp.h>

#include <algorithm>
#include <vector>

#include "zhurin_i_gauss_kernel/common/include/common.hpp"

namespace zhurin_i_gauss_kernel {

static int GetPixelMirror(const std::vector<std::vector<int>> &img, int row, int col, int width, int height) {
  if (row < 0) {
    row = -row - 1;
  } else if (row >= height) {
    row = 2 * height - row - 1;
  }
  if (col < 0) {
    col = -col - 1;
  } else if (col >= width) {
    col = 2 * width - col - 1;
  }
  return img[row][col];
}

ZhurinIGaussKernelALL::ZhurinIGaussKernelALL(const InType &in) {
  SetTypeOfTask(GetStaticTypeOfTask());
  GetInput() = in;
  GetOutput() = OutType{};
}

bool ZhurinIGaussKernelALL::ValidationImpl() {
  const auto &in = GetInput();
  int w = std::get<0>(in);
  int h = std::get<1>(in);
  int parts = std::get<2>(in);
  const auto &img = std::get<3>(in);

  if (w <= 0 || h <= 0 || parts <= 0 || parts > w) {
    return false;
  }
  if (std::cmp_not_equal(img.size(), h)) {
    return false;
  }
  for (int i = 0; i < h; ++i) {
    if (std::cmp_not_equal(img[i].size(), w)) {
      return false;
    }
  }

  int initialized = 0;
  MPI_Initialized(&initialized);
  if (!initialized) {
    return false;
  }

  return true;
}

bool ZhurinIGaussKernelALL::PreProcessingImpl() {
  const auto &in = GetInput();
  width_ = std::get<0>(in);
  height_ = std::get<1>(in);
  num_parts_ = std::get<2>(in);
  image_ = std::get<3>(in);
  result_.assign(height_, std::vector<int>(width_, 0));
  return true;
}

bool ZhurinIGaussKernelALL::RunImpl() {
  int rank, size;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &size);

  const int w = width_;
  const int h = height_;
  const auto &img = image_;

  int base = h / size;
  int rem = h % size;
  int start = rank * base + std::min(rank, rem);
  int end = start + base + (rank < rem ? 1 : 0);
  int local_rows = end - start;

  std::vector<int> flat_result(h * w, 0);
  std::vector<int> local_flat(local_rows * w, 0);

#pragma omp parallel for default(none) shared(local_flat, w, start, end, img, h)
  for (int i = start; i < end; ++i) {
    int row_idx = i - start;
    for (int j = 0; j < w; ++j) {
      int sum = 0;
      for (int ki = -1; ki <= 1; ++ki) {
        for (int kj = -1; kj <= 1; ++kj) {
          sum += GetPixelMirror(img, i + ki, j + kj, w, h) * kKernel[ki + 1][kj + 1];
        }
      }
      local_flat[row_idx * w + j] = sum >> kShift;
    }
  }

  std::vector<int> recv_counts(size, 0);
  std::vector<int> displs(size, 0);
  MPI_Allgather(&local_rows, 1, MPI_INT, recv_counts.data(), 1, MPI_INT, MPI_COMM_WORLD);
  int total = 0;
  for (int i = 0; i < size; ++i) {
    displs[i] = total;
    total += recv_counts[i];
  }
  MPI_Allgatherv(local_flat.data(), local_rows * w, MPI_INT, flat_result.data(), recv_counts.data(), displs.data(),
                 MPI_INT, MPI_COMM_WORLD);

  for (int i = 0; i < h; ++i) {
    std::copy(flat_result.begin() + i * w, flat_result.begin() + (i + 1) * w, result_[i].begin());
  }

  return true;
}

bool ZhurinIGaussKernelALL::PostProcessingImpl() {
  GetOutput() = std::move(result_);
  return true;
}

}  // namespace zhurin_i_gauss_kernel
