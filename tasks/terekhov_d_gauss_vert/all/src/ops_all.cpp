#include "terekhov_d_gauss_vert/all/include/ops_all.hpp"

#include <mpi.h>
#include <omp.h>

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <vector>

#include "terekhov_d_gauss_vert/common/include/common.hpp"

namespace terekhov_d_gauss_vert {

namespace {

inline int Mirror(int coord, int max_val) {
  if (coord < 0) {
    return -coord - 1;
  }
  if (coord >= max_val) {
    return (2 * max_val) - coord - 1;
  }
  return coord;
}

inline void FillLocalPadded(std::vector<int> &local_padded, const InType &in, int padded_width, int width, int height,
                            int start_row, int local_padded_height) {
  for (int row = 0; row < local_padded_height; ++row) {
    const int global_row = start_row + row - 1;
    for (int col = 0; col < padded_width; ++col) {
      const int src_x = Mirror(col - 1, width);
      const int src_y = Mirror(global_row, height);
      const size_t src_idx = (static_cast<size_t>(src_y) * static_cast<size_t>(width)) + static_cast<size_t>(src_x);
      const size_t padded_idx =
          (static_cast<size_t>(row) * static_cast<size_t>(padded_width)) + static_cast<size_t>(col);
      local_padded[padded_idx] = in.data[src_idx];
    }
  }
}

inline void ProcessLocalRows(std::vector<int> &local_result, const std::vector<int> &local_padded, int padded_width,
                             int width, int local_height) {
#pragma omp parallel for schedule(static) default(none) \
    shared(local_result, local_padded, padded_width, width, local_height) shared(kGaussKernel)
  for (int row = 0; row < local_height; ++row) {
    for (int col = 0; col < width; ++col) {
      const size_t idx = (static_cast<size_t>(row) * static_cast<size_t>(width)) + static_cast<size_t>(col);
      float sum = 0.0F;

      for (int ky = -1; ky <= 1; ++ky) {
        for (int kx = -1; kx <= 1; ++kx) {
          const int px = col + kx + 1;
          const int py = row + ky + 1;
          const int kernel_idx = ((ky + 1) * 3) + (kx + 1);
          const size_t padded_idx =
              (static_cast<size_t>(py) * static_cast<size_t>(padded_width)) + static_cast<size_t>(px);
          sum += static_cast<float>(local_padded[padded_idx]) * kGaussKernel[static_cast<size_t>(kernel_idx)];
        }
      }

      local_result[idx] = static_cast<int>(std::lround(sum));
    }
  }
}

inline void GatherResults(std::vector<int> &output_data, const std::vector<int> &local_result, int size) {
  std::vector<int> recv_counts(static_cast<size_t>(size));
  std::vector<int> displs(static_cast<size_t>(size));

  int local_count = static_cast<int>(local_result.size());
  MPI_Allgather(&local_count, 1, MPI_INT, recv_counts.data(), 1, MPI_INT, MPI_COMM_WORLD);

  displs[0] = 0;
  for (int idx = 1; idx < size; ++idx) {
    displs[idx] = displs[idx - 1] + recv_counts[idx - 1];
  }

  MPI_Allgatherv(local_result.data(), local_count, MPI_INT, output_data.data(), recv_counts.data(), displs.data(),
                 MPI_INT, MPI_COMM_WORLD);
}

inline OutType SolveALL(const InType &in) {
  int rank = 0;
  int size = 1;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &size);

  const int width = in.width;
  const int height = in.height;

  OutType output;
  output.width = width;
  output.height = height;
  output.data.resize(static_cast<size_t>(width) * static_cast<size_t>(height));

  int rows_per_proc = height / size;
  int remainder = height % size;
  int start_row = (rank * rows_per_proc) + std::min(rank, remainder);
  int end_row = start_row + rows_per_proc + ((rank < remainder) ? 1 : 0);

  if (start_row >= end_row) {
    return output;
  }

  const int padded_width = width + 2;
  const int local_height = end_row - start_row;
  const int local_padded_height = local_height + 2;

  std::vector<int> local_padded(static_cast<size_t>(padded_width) * static_cast<size_t>(local_padded_height), 0);
  FillLocalPadded(local_padded, in, padded_width, width, height, start_row, local_padded_height);

  std::vector<int> local_result(static_cast<size_t>(local_height) * static_cast<size_t>(width), 0);
  ProcessLocalRows(local_result, local_padded, padded_width, width, local_height);

  GatherResults(output.data, local_result, size);

  return output;
}

}  // namespace

TerekhovDGaussVertALL::TerekhovDGaussVertALL(const InType &in) {
  SetTypeOfTask(GetStaticTypeOfTask());
  GetInput() = in;
  local_out_.width = 0;
  local_out_.height = 0;
  local_out_.data.clear();
}

bool TerekhovDGaussVertALL::ValidationImpl() {
  const auto &in = GetInput();
  if ((in.width <= 0) || (in.height <= 0)) {
    return false;
  }
  const std::size_t need = static_cast<std::size_t>(in.width) * static_cast<std::size_t>(in.height);
  return in.data.size() == need;
}

bool TerekhovDGaussVertALL::PreProcessingImpl() {
  local_out_.width = 0;
  local_out_.height = 0;
  local_out_.data.clear();
  return true;
}

bool TerekhovDGaussVertALL::RunImpl() {
  if (!ValidationImpl()) {
    return false;
  }
  local_out_ = SolveALL(GetInput());
  return true;
}

bool TerekhovDGaussVertALL::PostProcessingImpl() {
  GetOutput() = local_out_;
  const auto &out = GetOutput();
  return out.data.size() == (static_cast<size_t>(out.width) * static_cast<size_t>(out.height));
}

}  // namespace terekhov_d_gauss_vert
