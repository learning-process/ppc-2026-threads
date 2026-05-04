#include "zhurin_i_gauss_kernel/all/include/ops_all.hpp"

#include <mpi.h>
#include <omp.h>

#include <algorithm>
#include <cstddef>
#include <utility>
#include <vector>

#include "zhurin_i_gauss_kernel/common/include/common.hpp"

namespace zhurin_i_gauss_kernel {

namespace {

int GetPixelMirror(const std::vector<std::vector<int>> &img, int row, int col, int width, int height) {
  if (row < 0) {
    row = -row - 1;
  } else if (row >= height) {
    row = (2 * height) - row - 1;
  }
  if (col < 0) {
    col = -col - 1;
  } else if (col >= width) {
    col = (2 * width) - col - 1;
  }
  return img[row][col];
}

}  // namespace

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
  if (std::cmp_not_equal(img.size(), static_cast<std::size_t>(h))) {
    return false;
  }
  for (int i = 0; i < h; ++i) {
    if (std::cmp_not_equal(img[i].size(), static_cast<std::size_t>(w))) {
      return false;
    }
  }

  int initialized = 0;
  MPI_Initialized(&initialized);
  return initialized != 0;
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
  int rank = 0;
  int size = 1;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &size);

  const int w = width_;
  const int h = height_;
  const auto &img = image_;

  const int base = h / size;
  const int rem = h % size;
  const int start = (rank * base) + std::min(rank, rem);
  const int end = start + base + ((rank < rem) ? 1 : 0);
  const int local_rows = end - start;

  std::vector<int> flat_result(static_cast<std::size_t>(h) * static_cast<std::size_t>(w), 0);
  std::vector<int> local_flat(static_cast<std::size_t>(local_rows) * static_cast<std::size_t>(w), 0);

#pragma omp parallel for default(none) shared(local_flat, w, start, end, img, h)
  for (int i = start; i < end; ++i) {
    const int row_idx = i - start;
    for (int j = 0; j < w; ++j) {
      const int p00 = GetPixelMirror(img, i - 1, j - 1, w, h);
      const int p01 = GetPixelMirror(img, i - 1, j, w, h);
      const int p02 = GetPixelMirror(img, i - 1, j + 1, w, h);
      const int p10 = GetPixelMirror(img, i, j - 1, w, h);
      const int p11 = GetPixelMirror(img, i, j, w, h);
      const int p12 = GetPixelMirror(img, i, j + 1, w, h);
      const int p20 = GetPixelMirror(img, i + 1, j - 1, w, h);
      const int p21 = GetPixelMirror(img, i + 1, j, w, h);
      const int p22 = GetPixelMirror(img, i + 1, j + 1, w, h);

      int sum = (p00 * kKernel[0][0]) + (p01 * kKernel[0][1]) + (p02 * kKernel[0][2]) + (p10 * kKernel[1][0]) +
                (p11 * kKernel[1][1]) + (p12 * kKernel[1][2]) + (p20 * kKernel[2][0]) + (p21 * kKernel[2][1]) +
                (p22 * kKernel[2][2]);

      const std::size_t idx =
          static_cast<std::size_t>(row_idx) * static_cast<std::size_t>(w) + static_cast<std::size_t>(j);
      local_flat[idx] = sum >> kShift;
    }
  }

  std::vector<int> recv_counts(static_cast<std::size_t>(size), 0);
  std::vector<int> displs(static_cast<std::size_t>(size), 0);
  MPI_Allgather(&local_rows, 1, MPI_INT, recv_counts.data(), 1, MPI_INT, MPI_COMM_WORLD);
  int total = 0;
  for (int i = 0; i < size; ++i) {
    displs[i] = total;
    total += recv_counts[i];
  }
  MPI_Allgatherv(local_flat.data(), local_rows * w, MPI_INT, flat_result.data(), recv_counts.data(), displs.data(),
                 MPI_INT, MPI_COMM_WORLD);

  for (int i = 0; i < h; ++i) {
    const std::ptrdiff_t offset = static_cast<std::ptrdiff_t>(i) * static_cast<std::ptrdiff_t>(w);
    std::copy_n(flat_result.data() + offset, static_cast<std::ptrdiff_t>(w), result_[i].begin());
  }

  return true;
}

bool ZhurinIGaussKernelALL::PostProcessingImpl() {
  GetOutput() = std::move(result_);
  return true;
}

}  // namespace zhurin_i_gauss_kernel
