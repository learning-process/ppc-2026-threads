#include "boltenkov_s_gaussian_kernel/all/include/ops_all.hpp"

#include <mpi.h>
#include <omp.h>

#include <algorithm>
#include <climits>
#include <cstddef>
#include <vector>

#include "boltenkov_s_gaussian_kernel/common/include/common.hpp"
#include "util/include/util.hpp"

namespace boltenkov_s_gaussian_kernel {

BoltenkovSGaussianKernelALL::BoltenkovSGaussianKernelALL(const InType &in) {
  SetTypeOfTask(GetStaticTypeOfTask());
  int rank = 0;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  if (rank == 0) {
    GetInput() = in;
  } else {
    GetInput() = InType();
  }
  GetOutput() = std::vector<std::vector<int>>();
}

bool BoltenkovSGaussianKernelALL::ValidationImpl() {
  int rank = 0;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  if (rank == 0) {
    std::size_t n = std::get<0>(GetInput());
    std::size_t m = std::get<1>(GetInput());
    const auto &data = std::get<2>(GetInput());
    if (data.size() != n) {
      return false;
    }
    for (std::size_t i = 0; i < n; ++i) {
      if (data[i].size() != m) {
        return false;
      }
    }
    return true;
  }
  return true;
}

bool BoltenkovSGaussianKernelALL::PreProcessingImpl() {
  int rank = 0;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  auto n_size_t = std::get<0>(GetInput());
  auto m_size_t = std::get<1>(GetInput());
  if (n_size_t > INT_MAX || m_size_t > INT_MAX) {
    return false;
  }
  int n_val = 0;
  int m_val = 0;
  if (rank == 0) {
    n_val = static_cast<int>(n_size_t);
    m_val = static_cast<int>(m_size_t);
  }
  MPI_Bcast(&n_val, 1, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Bcast(&m_val, 1, MPI_INT, 0, MPI_COMM_WORLD);
  if (n_val < 1e6 && n_val > 0) {
    GetOutput().resize(static_cast<std::size_t>(n_val));
  } else {
    return false;
  }
  for (int i = 0; i < n_val; ++i) {
    if (m_val < 1e6 && m_val > 0) {
      GetOutput()[i].resize(static_cast<std::size_t>(m_val));
    } else {
      return false;
    }
  }
  return true;
}

bool BoltenkovSGaussianKernelALL::IsValidSize(int n, int m) {
  return n > 0 && m > 0 && n < 1e6 && m < 1e6;
}

void BoltenkovSGaussianKernelALL::ComputeScatterParams(int n, int m, int size, int rows_per_proc,
                                                       std::vector<int> &send_counts, std::vector<int> &displs) {
  send_counts.assign(size, 0);
  displs.assign(size, 0);
  for (int i = 0; i < size; ++i) {
    int s = i * rows_per_proc;
    int e = std::min(s + rows_per_proc, n) - 1;
    int rows = (s < n) ? (e - s + 1) : 0;
    if (rows > 0) {
      int h_first = std::max(0, s - 1);
      int h_last = std::min(n - 1, e + 1);
      int h_rows = h_last - h_first + 1;
      send_counts[i] = h_rows * m;
      displs[i] = h_first * m;
    }
  }
}

void BoltenkovSGaussianKernelALL::ComputeGatherDispls(int m, const std::vector<int> &gather_counts,
                                                      std::vector<int> &recv_counts, std::vector<int> &recv_displs) {
  int size = static_cast<int>(gather_counts.size());
  recv_counts.resize(size);
  recv_displs.resize(size);
  for (int i = 0; i < size; ++i) {
    recv_counts[i] = gather_counts[i] * m;
    recv_displs[i] = (i == 0) ? 0 : recv_displs[i - 1] + recv_counts[i - 1];
  }
}

std::vector<int> BoltenkovSGaussianKernelALL::ApplyGaussianFilterFlat(const std::vector<int> &local_halo_flat,
                                                                      int halo_rows, int local_start_row,
                                                                      int local_rows, int m, const int kernel[3][3],
                                                                      int shift) {
  const int tmp_rows = local_rows + 2;
  const int tmp_cols = m + 2;
  std::vector<int> tmp(static_cast<size_t>(tmp_rows) * static_cast<size_t>(tmp_cols), 0);

  const int halo_first = std::max(0, local_start_row - 1);

  for (int i = 0; i < tmp_rows; ++i) {
    int global_row = local_start_row - 1 + i;
    if (global_row >= halo_first && global_row < halo_first + halo_rows) {
      const int src_offset = (global_row - halo_first) * m;
      int *dst_row = &tmp[static_cast<size_t>(i) * static_cast<size_t>(tmp_cols) + 1];
      std::copy_n(&local_halo_flat[src_offset], m, dst_row);
    }
  }

  std::vector<int> local_res(static_cast<size_t>(local_rows) * static_cast<size_t>(m), 0);

#pragma omp parallel for num_threads(ppc::util::GetNumThreads()) default(none) \
    shared(tmp, local_res, local_rows, m, kernel, shift, tmp_cols)
  for (int i = 0; i < local_rows; ++i) {
    const int *row0 = &tmp[static_cast<size_t>(i) * static_cast<size_t>(tmp_cols)];
    const int *row1 = row0 + tmp_cols;
    const int *row2 = row1 + tmp_cols;
    int *out_row = &local_res[static_cast<size_t>(i) * static_cast<size_t>(m)];

    const int k00 = kernel[0][0];
    const int k01 = kernel[0][1];
    const int k02 = kernel[0][2];
    const int k10 = kernel[1][0];
    const int k11 = kernel[1][1];
    const int k12 = kernel[1][2];
    const int k20 = kernel[2][0];
    const int k21 = kernel[2][1];
    const int k22 = kernel[2][2];

    for (int j = 0; j < m; ++j) {
      int val = (row0[j] * k00) + (row0[j + 1] * k01) + (row0[j + 2] * k02) + (row1[j] * k10) + (row1[j + 1] * k11) +
                (row1[j + 2] * k12) + (row2[j] * k20) + (row2[j + 1] * k21) + (row2[j + 2] * k22);
      out_row[j] = val >> shift;
    }
  }

  return local_res;
}

bool BoltenkovSGaussianKernelALL::RunImpl() {
  int rank = 0;
  int size = 0;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &size);

  int n = static_cast<int>(GetOutput().size());
  int m = static_cast<int>(GetOutput()[0].size());
  if (!IsValidSize(n, m)) {
    return false;
  }

  MPI_Bcast(&n, 1, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Bcast(&m, 1, MPI_INT, 0, MPI_COMM_WORLD);

  std::vector<int> data_flat(static_cast<size_t>(n) * static_cast<size_t>(m));
  if (rank == 0) {
    const auto &global_data = std::get<2>(GetInput());
    for (int i = 0; i < n; ++i) {
      std::copy_n(global_data[i].data(), m, &data_flat[static_cast<size_t>(i) * static_cast<size_t>(m)]);
    }
  }

  int rows_per_proc = (n + size - 1) / size;
  int local_start = rank * rows_per_proc;
  int local_rows = 0;
  if (local_start < n) {
    local_rows = std::min(rows_per_proc, n - local_start);
  }

  int halo_first = std::max(0, local_start - 1);
  int halo_last = std::min(n - 1, local_start + local_rows);
  int halo_rows = (local_rows > 0) ? (halo_last - halo_first + 1) : 0;

  std::vector<int> send_counts(size, 0);
  std::vector<int> displs(size, 0);
  if (rank == 0) {
    ComputeScatterParams(n, m, size, rows_per_proc, send_counts, displs);
  }

  std::vector<int> local_halo_flat(static_cast<size_t>(halo_rows) * static_cast<size_t>(m));

  MPI_Scatterv(data_flat.data(), send_counts.data(), displs.data(), MPI_INT, local_halo_flat.data(),
               static_cast<int>(local_halo_flat.size()), MPI_INT, 0, MPI_COMM_WORLD);

  std::vector<int> local_res_flat;
  if (local_rows > 0) {
    local_res_flat = ApplyGaussianFilterFlat(local_halo_flat, halo_rows, local_start, local_rows, m, kernel_, shift_);
  }

  std::vector<int> gather_counts(size, 0);
  MPI_Gather(&local_rows, 1, MPI_INT, gather_counts.data(), 1, MPI_INT, 0, MPI_COMM_WORLD);

  std::vector<int> recv_counts(size, 0);
  std::vector<int> recv_displs(size, 0);
  if (rank == 0) {
    ComputeGatherDispls(m, gather_counts, recv_counts, recv_displs);
  }

  std::vector<int> out_flat(static_cast<size_t>(n) * static_cast<size_t>(m));

  MPI_Gatherv(local_res_flat.data(), static_cast<int>(local_res_flat.size()), MPI_INT, out_flat.data(),
              recv_counts.data(), recv_displs.data(), MPI_INT, 0, MPI_COMM_WORLD);

  MPI_Bcast(out_flat.data(), static_cast<int>(out_flat.size()), MPI_INT, 0, MPI_COMM_WORLD);

  auto &output = GetOutput();
  for (int i = 0; i < n; ++i) {
    std::copy_n(&out_flat[static_cast<size_t>(i) * static_cast<size_t>(m)], m, output[i].data());
  }

  return true;
}

bool BoltenkovSGaussianKernelALL::PostProcessingImpl() {
  return true;
}

}  // namespace boltenkov_s_gaussian_kernel
