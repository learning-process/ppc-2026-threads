#include "romanova_v_linear_histogram_stretch/all/include/ops_all.hpp"

#include <mpi.h>

#include <algorithm>
#include <cstddef>
#include <cstdint>
#include <utility>
#include <vector>

#include "romanova_v_linear_histogram_stretch/common/include/common.hpp"

namespace romanova_v_linear_histogram_stretch_threads {

RomanovaVLinHistogramStretchALL::RomanovaVLinHistogramStretchALL(const InType &in) {
  SetTypeOfTask(GetStaticTypeOfTask());
  GetInput() = in;
  GetOutput().clear();
}

bool RomanovaVLinHistogramStretchALL::ValidationImpl() {
  bool status = true;
  int rank = 0;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  if (rank == 0) {
    status = !GetInput().empty();
  }

  MPI_Bcast(&status, 1, MPI_C_BOOL, 0, MPI_COMM_WORLD);
  return status;
}

bool RomanovaVLinHistogramStretchALL::PreProcessingImpl() {
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

  size_t delta = size / n;
  size_t extra = size % n;

  local_size_ = delta + (std::cmp_less(rank, extra) ? 1 : 0);

  vector_counts_.resize(n);
  vector_displs_.resize(n);

  if (rank == 0) {
    for (int i = 0; i < n; i++) {
      vector_counts_[i] = static_cast<int>(delta + (std::cmp_less(i, extra) ? 1 : 0));
    }

    for (int i = 1; i < n; i++) {
      vector_displs_[i] = vector_displs_[i - 1] + vector_counts_[i - 1];
    }
  }

  MPI_Bcast(vector_counts_.data(), n, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Bcast(vector_displs_.data(), n, MPI_INT, 0, MPI_COMM_WORLD);

  local_data_.resize(local_size_);

  MPI_Scatterv(rank == 0 ? GetInput().data() : nullptr, vector_counts_.data(), vector_displs_.data(), MPI_UNSIGNED_CHAR,
               local_data_.data(), static_cast<int>(local_size_), MPI_UNSIGNED_CHAR, 0, MPI_COMM_WORLD);

  local_out_ = OutType(local_size_);
  return !GetOutput().empty();
}

bool RomanovaVLinHistogramStretchALL::RunImpl() {
  uint8_t loc_min_v = 255;
  uint8_t loc_max_v = 0;

#pragma omp parallel for default(none) /*shared(local_data_, local_size_)*/ reduction(min : loc_min_v) \
    reduction(max : loc_max_v)
  for (size_t i = 0; i < local_size_; i++) {
    loc_min_v = std::min(loc_min_v, local_data_[i]);
    loc_max_v = std::max(loc_max_v, local_data_[i]);
  }

  uint8_t min_v = 255;
  uint8_t max_v = 0;

  MPI_Allreduce(static_cast<void *>(&loc_min_v), static_cast<void *>(&min_v), 1, MPI_UNSIGNED_CHAR, MPI_MIN,
                MPI_COMM_WORLD);
  MPI_Allreduce(static_cast<void *>(&loc_max_v), static_cast<void *>(&max_v), 1, MPI_UNSIGNED_CHAR, MPI_MAX,
                MPI_COMM_WORLD);

  if (min_v == max_v) {
#pragma omp parallel for default(none)  // shared(local_out_, local_data_, local_size_)
    for (size_t i = 0; i < local_size_; i++) {
      local_out_[i] = local_data_[i];
    }
    return true;
  }

  const uint8_t diff = max_v - min_v;
  const double ratio = 255.0 / diff;

#pragma omp parallel for default(none) shared(min_v, ratio)
  for (size_t i = 0; i < local_size_; i++) {
    uint8_t pix = local_data_[i];
    local_out_[i] = (std::clamp(static_cast<int>((pix - min_v) * ratio), 0, 255));
  }

  return true;
}

bool RomanovaVLinHistogramStretchALL::PostProcessingImpl() {
  MPI_Allgatherv(local_out_.data(), static_cast<int>(local_size_), MPI_UNSIGNED_CHAR, GetOutput().data(),
                 vector_counts_.data(), vector_displs_.data(), MPI_UNSIGNED_CHAR, MPI_COMM_WORLD);
  return true;
}

}  // namespace romanova_v_linear_histogram_stretch_threads
