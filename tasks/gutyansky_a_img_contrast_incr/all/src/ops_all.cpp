#include "gutyansky_a_img_contrast_incr/all/include/ops_all.hpp"

#include <mpi.h>

#include <algorithm>
#include <array>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <utility>
#include <vector>

#include "gutyansky_a_img_contrast_incr/common/include/common.hpp"

namespace gutyansky_a_img_contrast_incr {

GutyanskyAImgContrastIncrALL::GutyanskyAImgContrastIncrALL(const InType &in) {
  SetTypeOfTask(GetStaticTypeOfTask());
  GetInput() = in;
  GetOutput() = {};
}

bool GutyanskyAImgContrastIncrALL::ValidationImpl() {
  return !GetInput().empty();
}

bool GutyanskyAImgContrastIncrALL::PreProcessingImpl() {
  return true;
}

bool GutyanskyAImgContrastIncrALL::RunImpl() {
  int rank = -1;
  int p_count = -1;

  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &p_count);

  size_t sz = GetInput().size();

  MPI_Bcast(&sz, 1, MPI_UINT64_T, 0, MPI_COMM_WORLD);

  size_t chunk_size = sz / static_cast<size_t>(p_count);
  size_t remainder_size = sz % static_cast<size_t>(p_count);

  size_t local_count = chunk_size + (std::cmp_less(rank, remainder_size) ? 1 : 0);

  std::vector<uint8_t> local_data(local_count);

  std::vector<int> send_counts(p_count);
  std::vector<int> locs(p_count);

  for (int i = 0; i < p_count; i++) {
    size_t elems_for_proc = chunk_size + (std::cmp_less(i, remainder_size) ? 1 : 0);
    send_counts[i] = static_cast<int>(elems_for_proc);
  }

  for (int i = 1; i < p_count; i++) {
    locs[i] = locs[i - 1] + send_counts[i - 1];
  }

  MPI_Scatterv(GetInput().data(), send_counts.data(), locs.data(), MPI_UINT8_T, local_data.data(),
               static_cast<int>(local_count), MPI_UINT8_T, 0, MPI_COMM_WORLD);

  uint8_t local_lower_bound = std::numeric_limits<uint8_t>::max();
  uint8_t local_upper_bound = std::numeric_limits<uint8_t>::min();

#pragma omp parallel for default(none) shared(local_data, local_count) reduction(min : local_lower_bound) \
    reduction(max : local_upper_bound)
  for (size_t i = 0; i < local_count; i++) {
    uint8_t val = local_data[i];
    local_lower_bound = std::min(local_lower_bound, val);
    local_upper_bound = std::max(local_upper_bound, val);
  }

  uint8_t lower_bound = std::numeric_limits<uint8_t>::max();
  uint8_t upper_bound = std::numeric_limits<uint8_t>::min();

  MPI_Allreduce(&local_lower_bound, &lower_bound, 1, MPI_UINT8_T, MPI_MIN, MPI_COMM_WORLD);
  MPI_Allreduce(&local_upper_bound, &upper_bound, 1, MPI_UINT8_T, MPI_MAX, MPI_COMM_WORLD);

  uint8_t delta = upper_bound - lower_bound;

  if (delta != 0) {
    std::array<uint8_t, 256> lut{};
    for (int pix = 0; pix < 256; pix++) {
      int shifted = pix - lower_bound;
      lut.at(pix) = static_cast<uint8_t>(255 * shifted / delta);
    }

#pragma omp parallel for default(none) shared(local_data, local_count, lut)
    for (size_t i = 0; i < local_count; i++) {
      local_data[i] = lut.at(local_data[i]);
    }
  }

  GetOutput().resize(sz);
  MPI_Allgatherv(local_data.data(), static_cast<int>(local_count), MPI_UINT8_T, GetOutput().data(), send_counts.data(),
                 locs.data(), MPI_UINT8_T, MPI_COMM_WORLD);

  return true;
}

bool GutyanskyAImgContrastIncrALL::PostProcessingImpl() {
  return true;
}

}  // namespace gutyansky_a_img_contrast_incr
