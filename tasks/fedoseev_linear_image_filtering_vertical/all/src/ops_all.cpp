#include "fedoseev_linear_image_filtering_vertical/all/include/ops_all.hpp"

#include <mpi.h>
#include <omp.h>

#include <algorithm>
#include <array>
#include <cstddef>
#include <vector>

#include "fedoseev_linear_image_filtering_vertical/common/include/common.hpp"
#include "util/include/util.hpp"

namespace fedoseev_linear_image_filtering_vertical {

namespace {
const std::array<std::array<int, 3>, 3> kKernel = {{{{1, 2, 1}}, {{2, 4, 2}}, {{1, 2, 1}}}};
const int kKernelSum = 16;

void BroadcastAndFilter(int w, int h, const std::vector<int> &img, std::vector<int> &result, int num_threads) {
#pragma omp parallel for num_threads(num_threads) default(none) shared(img, result, w, h, kKernel, kKernelSum)
  for (int row = 0; row < h; ++row) {
    for (int col = 0; col < w; ++col) {
      int sum = 0;
      for (int ky = -1; ky <= 1; ++ky) {
        for (int kx = -1; kx <= 1; ++kx) {
          int px = col + kx;
          int py = row + ky;
          px = std::clamp(px, 0, w - 1);
          py = std::clamp(py, 0, h - 1);
          sum += img[py * w + px] * kKernel[ky + 1][kx + 1];
        }
      }
      result[row * w + col] = sum / kKernelSum;
    }
  }
}

void ExchangeGhostRows(int rank, int size, int w, int local_rows, const std::vector<int> &local_data,
                       std::vector<int> &ghost_src) {
  bool is_mpi = ppc::util::IsUnderMpirun();
  if (!is_mpi) {
    std::copy(local_data.begin(), local_data.begin() + w, ghost_src.begin());
    std::copy(local_data.begin() + (local_rows - 1) * w, local_data.begin() + local_rows * w,
              ghost_src.begin() + (local_rows + 1) * w);
    return;
  }

  if (rank > 0) {
    MPI_Sendrecv(const_cast<int *>(local_data.data()), w, MPI_INT, rank - 1, 0, ghost_src.data(), w, MPI_INT, rank - 1,
                 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
  } else {
    std::copy(local_data.begin(), local_data.begin() + w, ghost_src.begin());
  }

  if (rank < size - 1) {
    MPI_Sendrecv(const_cast<int *>(local_data.data() + (local_rows - 1) * w), w, MPI_INT, rank + 1, 0,
                 ghost_src.data() + (local_rows + 1) * w, w, MPI_INT, rank + 1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
  } else {
    std::copy(local_data.begin() + (local_rows - 1) * w, local_data.begin() + local_rows * w,
              ghost_src.begin() + (local_rows + 1) * w);
  }
}
}  // namespace

LinearImageFilteringVerticalAll::LinearImageFilteringVerticalAll(const InType &in) {
  SetTypeOfTask(GetStaticTypeOfTask());
  GetInput() = in;
  GetOutput() = InType{};
}

bool LinearImageFilteringVerticalAll::ValidationImpl() {
  const InType &input = GetInput();
  if (input.width < 3 || input.height < 3) {
    return false;
  }
  return input.data.size() == static_cast<size_t>(input.width) * static_cast<size_t>(input.height);
}

bool LinearImageFilteringVerticalAll::PreProcessingImpl() {
  const InType &input = GetInput();
  OutType output;
  output.width = input.width;
  output.height = input.height;
  output.data.resize(static_cast<size_t>(input.width) * static_cast<size_t>(input.height), 0);
  GetOutput() = output;
  return true;
}

void LinearImageFilteringVerticalAll::DistributeData(int rank, int size, int &h) {
  bool is_mpi = ppc::util::IsUnderMpirun();
  if (is_mpi) {
    MPI_Bcast(&h, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(&w_, 1, MPI_INT, 0, MPI_COMM_WORLD);
  }
  if (h == 0) {
    return;
  }

  counts_.resize(size);
  displs_.resize(size);
  int rows_per_proc = h / size;
  int rem = h % size;
  int curr = 0;
  for (int i = 0; i < size; ++i) {
    int rows = rows_per_proc + (i < rem ? 1 : 0);
    counts_[i] = rows * w_;
    displs_[i] = curr;
    curr += counts_[i];
  }
  local_data_.resize(static_cast<size_t>(counts_[rank]));

  if (is_mpi) {
    MPI_Scatterv(rank == 0 ? GetInput().data.data() : nullptr, counts_.data(), displs_.data(), MPI_INT,
                 local_data_.data(), counts_[rank], MPI_INT, 0, MPI_COMM_WORLD);
  } else {
    local_data_ = GetInput().data;
  }
}

void LinearImageFilteringVerticalAll::LocalProcessing(int rank, int size, int num_threads) {
  int local_rows = static_cast<int>(local_data_.size()) / w_;
  if (local_rows == 0) {
    return;
  }

  int ghost_rows = local_rows + 2;
  std::vector<int> ghost_src(static_cast<size_t>(ghost_rows) * static_cast<size_t>(w_), 0);
  for (int i = 0; i < local_rows; ++i) {
    std::copy(local_data_.begin() + static_cast<ptrdiff_t>(i) * w_,
              local_data_.begin() + static_cast<ptrdiff_t>(i + 1) * w_,
              ghost_src.begin() + static_cast<ptrdiff_t>(i + 1) * w_);
  }

  ExchangeGhostRows(rank, size, w_, local_rows, local_data_, ghost_src);

  local_result_.resize(static_cast<size_t>(local_rows) * static_cast<size_t>(w_), 0);
#pragma omp parallel for num_threads(num_threads) default(none) \
    shared(ghost_src, local_result_, local_rows, ghost_rows, w_, kKernel, kKernelSum)
  for (int i = 0; i < local_rows; ++i) {
    int row_in_ghost = i + 1;
    for (int col = 0; col < w_; ++col) {
      int sum = 0;
      for (int ky = -1; ky <= 1; ++ky) {
        for (int kx = -1; kx <= 1; ++kx) {
          int px = col + kx;
          int py = row_in_ghost + ky;
          if (px < 0) {
            px = 0;
          }
          if (px >= w_) {
            px = w_ - 1;
          }
          if (py < 0) {
            py = 0;
          }
          if (py >= ghost_rows) {
            py = ghost_rows - 1;
          }
          sum += ghost_src[py * w_ + px] * kKernel.at(ky + 1).at(kx + 1);
        }
      }
      local_result_[i * w_ + col] = sum / kKernelSum;
    }
  }
}

void LinearImageFilteringVerticalAll::GatherData(int rank, int /*size*/) {
  bool is_mpi = ppc::util::IsUnderMpirun();
  if (is_mpi) {
    MPI_Gatherv(local_result_.data(), static_cast<int>(local_result_.size()), MPI_INT,
                rank == 0 ? GetOutput().data.data() : nullptr, counts_.data(), displs_.data(), MPI_INT, 0,
                MPI_COMM_WORLD);
  } else {
    OutType out;
    out.width = w_;
    out.height = static_cast<int>(local_result_.size()) / w_;
    out.data = std::move(local_result_);
    GetOutput() = out;
  }
}

bool LinearImageFilteringVerticalAll::RunImpl() {
  int rank = 0, size = 1;
  bool is_mpi = ppc::util::IsUnderMpirun();
  if (is_mpi) {
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
  }

  const auto &input = GetInput();
  int w = input.width;
  int h = input.height;
  w_ = w;

  if (w * h <= 256) {
    std::vector<int> img;
    if (rank == 0) {
      img = input.data;
    }
    if (is_mpi) {
      MPI_Bcast(&w, 1, MPI_INT, 0, MPI_COMM_WORLD);
      MPI_Bcast(&h, 1, MPI_INT, 0, MPI_COMM_WORLD);
      if (rank != 0) {
        img.resize(w * h);
      }
      MPI_Bcast(img.data(), w * h, MPI_INT, 0, MPI_COMM_WORLD);
    } else {
      img = input.data;
    }
    std::vector<int> result(w * h, 0);
    int num_threads = ppc::util::GetNumThreads();
    BroadcastAndFilter(w, h, img, result, num_threads);
    GetOutput().width = w;
    GetOutput().height = h;
    GetOutput().data = std::move(result);
    return true;
  }

  DistributeData(rank, size, h);
  if (h == 0) {
    return true;
  }
  int num_threads = ppc::util::GetNumThreads();
  LocalProcessing(rank, size, num_threads);
  GatherData(rank, size);
  return true;
}

bool LinearImageFilteringVerticalAll::PostProcessingImpl() {
  return true;
}

}  // namespace fedoseev_linear_image_filtering_vertical
