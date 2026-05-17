#include "pylaeva_s_inc_contrast_img_by_lsh/all/include/ops_all.hpp"

#include <mpi.h>

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <utility>
#include <vector>

#include "pylaeva_s_inc_contrast_img_by_lsh/common/include/common.hpp"

namespace pylaeva_s_inc_contrast_img_by_lsh {

PylaevaSIncContrastImgByLshALL::PylaevaSIncContrastImgByLshALL(const InType &in) {
  SetTypeOfTask(GetStaticTypeOfTask());
  GetInput() = in;
  GetOutput().clear();
}

bool PylaevaSIncContrastImgByLshALL::ValidationImpl() {
  bool status = true;
  int rank = 0;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  if (rank == 0) {
    status = !GetInput().empty();
  }

  MPI_Bcast(&status, 1, MPI_C_BOOL, 0, MPI_COMM_WORLD);
  return status;
}

bool PylaevaSIncContrastImgByLshALL::PreProcessingImpl() {
  int rank = 0;
  int n = 0;

  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &n);

  size_t size = 0;
  if (rank == 0) {
    size = GetInput().size();
  }
  MPI_Bcast(&size, 1, MPI_UNSIGNED_LONG_LONG, 0, MPI_COMM_WORLD);
  GetOutput().resize(size);

  size_t base_chunk = size / n;
  size_t remainder = size % n;

  local_size_ = base_chunk + (std::cmp_less(rank, remainder) ? 1 : 0);

  recv_counts_.resize(n);
  recv_displs_.resize(n);

  if (rank == 0) {
    for (int i = 0; i < n; i++) {
      recv_counts_[i] = static_cast<int>(base_chunk + (std::cmp_less(i, remainder) ? 1 : 0));
    }

    for (int i = 1; i < n; i++) {
      recv_displs_[i] = recv_displs_[i - 1] + recv_counts_[i - 1];
    }
  }

  MPI_Bcast(recv_counts_.data(), n, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Bcast(recv_displs_.data(), n, MPI_INT, 0, MPI_COMM_WORLD);

  local_data_.resize(local_size_);

  MPI_Scatterv(rank == 0 ? GetInput().data() : nullptr, recv_counts_.data(), recv_displs_.data(), MPI_UNSIGNED_CHAR,
               local_data_.data(), static_cast<int>(local_size_), MPI_UNSIGNED_CHAR, 0, MPI_COMM_WORLD);

  local_out_ = OutType(local_size_);
  return !GetOutput().empty();
}

bool PylaevaSIncContrastImgByLshALL::RunImpl() {
  uint8_t loc_min_v = 255;
  uint8_t loc_max_v = 0;

#pragma omp parallel for default(none) reduction(min : loc_min_v) reduction(max : loc_max_v)
  for (size_t i = 0; i < local_size_; i++) {
    loc_min_v = std::min(loc_min_v, local_data_[i]);
    loc_max_v = std::max(loc_max_v, local_data_[i]);
  }

  uint8_t global_min = 255;
  uint8_t global_max = 0;

  MPI_Allreduce(static_cast<void *>(&loc_min_v), static_cast<void *>(&global_min), 1, MPI_UNSIGNED_CHAR, MPI_MIN,
                MPI_COMM_WORLD);
  MPI_Allreduce(static_cast<void *>(&loc_max_v), static_cast<void *>(&global_max), 1, MPI_UNSIGNED_CHAR, MPI_MAX,
                MPI_COMM_WORLD);

  if (global_min == global_max) {
#pragma omp parallel for default(none)
    for (size_t i = 0; i < local_size_; i++) {
      local_out_[i] = local_data_[i];
    }
    return true;
  }

  const int diff = global_max - global_min;
  const double scale = 255.0 / diff;

#pragma omp parallel for default(none) shared(global_min, scale)
  for (size_t i = 0; i < local_size_; i++) {
    int pixel = local_data_[i];
    int result = static_cast<int>(std::lround((pixel - global_min) * scale));
    local_out_[i] = static_cast<uint8_t>(std::clamp(result, 0, 255));
  }

  return true;
}

bool PylaevaSIncContrastImgByLshALL::PostProcessingImpl() {
  MPI_Allgatherv(local_out_.data(), static_cast<int>(local_size_), MPI_UNSIGNED_CHAR, GetOutput().data(),
                 recv_counts_.data(), recv_displs_.data(), MPI_UNSIGNED_CHAR, MPI_COMM_WORLD);
  return true;
}

}  // namespace pylaeva_s_inc_contrast_img_by_lsh
