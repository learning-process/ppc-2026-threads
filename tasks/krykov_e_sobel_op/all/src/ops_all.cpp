#include "krykov_e_sobel_op/all/include/ops_all.hpp"

#include <mpi.h>
#include <omp.h>

#include <array>
#include <cmath>
#include <cstddef>
#include <thread>
#include <vector>

#include "krykov_e_sobel_op/common/include/common.hpp"

namespace krykov_e_sobel_op {

KrykovESobelOpALL::KrykovESobelOpALL(const InType &in) {
  SetTypeOfTask(GetStaticTypeOfTask());
  GetInput() = in;
  GetOutput().clear();
}

bool KrykovESobelOpALL::ValidationImpl() {
  const auto &img = GetInput();
  return img.width > 2 && img.height > 2 && static_cast<int>(img.data.size()) == img.width * img.height;
}

bool KrykovESobelOpALL::PreProcessingImpl() {
  const auto &img = GetInput();

  width_ = img.width;
  height_ = img.height;

  grayscale_.resize(static_cast<size_t>(width_) * static_cast<size_t>(height_));
  // RGB → grayscale
  for (int i = 0; i < width_ * height_; ++i) {
    const Pixel &p = img.data[i];
    grayscale_[i] = static_cast<int>((0.299 * p.r) + (0.587 * p.g) + (0.114 * p.b));
  }
  GetOutput().assign(static_cast<size_t>(width_) * static_cast<size_t>(height_), 0);
  return true;
}

bool KrykovESobelOpALL::RunImpl() {
  const std::array<std::array<int, 3>, 3> gx_kernel = {{{-1, 0, 1}, {-2, 0, 2}, {-1, 0, 1}}};
  const std::array<std::array<int, 3>, 3> gy_kernel = {{{-1, -2, -1}, {0, 0, 0}, {1, 2, 1}}};

  auto &output = GetOutput();
  const auto &gray = grayscale_;
  const int h = height_;
  const int w = width_;

  int rank = 0;
  int nproc = 1;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &nproc);

  const int base = (h - 2) / nproc;
  const int extra = (h - 2) % nproc;

  const int local_count = base + (rank < extra ? 1 : 0);
  const int local_start = 1 + base * rank + std::min(rank, extra);

  std::vector<int> local_output(static_cast<size_t>(local_count) * static_cast<size_t>(w), 0);

  std::vector<int> counts(nproc);
  std::vector<int> displs(nproc);
  for (int i = 0; i < nproc; ++i) {
    const int cnt = base + (i < extra ? 1 : 0);
    const int start = 1 + base * i + std::min(i, extra);
    counts[i] = cnt * w;
    displs[i] = start * w;
  }

#pragma omp parallel for default(none) shared(local_output, gray, gx_kernel, gy_kernel) \
    firstprivate(local_start, local_count, w) schedule(static)
  for (int li = 0; li < local_count; ++li) {
    const int row = local_start + li;
    for (int col = 1; col < w - 1; ++col) {
      int gx = 0;
      int gy = 0;

      for (int ky = -1; ky <= 1; ++ky) {
        for (int kx = -1; kx <= 1; ++kx) {
          int pixel = gray[static_cast<size_t>((row + ky) * w + (col + kx))];
          gx += pixel * gx_kernel.at(static_cast<size_t>(ky + 1)).at(static_cast<size_t>(kx + 1));
          gy += pixel * gy_kernel.at(static_cast<size_t>(ky + 1)).at(static_cast<size_t>(kx + 1));
        }
      }

      int magnitude = static_cast<int>(std::sqrt(static_cast<double>((gx * gx) + (gy * gy))));
      local_output[static_cast<size_t>(li * w + col)] = magnitude;
    }
  }

  MPI_Gatherv(local_output.data(), counts[rank], MPI_INT, output.data(), counts.data(), displs.data(), MPI_INT, 0,
              MPI_COMM_WORLD);
  MPI_Bcast(output.data(), w * h, MPI_INT, 0, MPI_COMM_WORLD);

  return true;
}

bool KrykovESobelOpALL::PostProcessingImpl() {
  return true;
}

}  // namespace krykov_e_sobel_op
