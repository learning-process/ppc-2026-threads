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

BoltenkovSGaussianKernelALL::BoltenkovSGaussianKernelALL(const InType &in)
    : kernel_{{{1, 2, 1}, {2, 4, 2}, {1, 2, 1}}} {
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

void BoltenkovSGaussianKernelALL::BcastSizes(int &n, int &m, int rank) {
  if (rank == 0) {
    n = static_cast<int>(std::get<0>(GetInput()));
    m = static_cast<int>(std::get<1>(GetInput()));
  }
  MPI_Bcast(&n, 1, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Bcast(&m, 1, MPI_INT, 0, MPI_COMM_WORLD);
}

void BoltenkovSGaussianKernelALL::SendRowsToOneProcess(int proc, int rows_per_proc, int n, int m,
                                                       const std::vector<std::vector<int>> &global_data) {
  int p_start = proc * rows_per_proc;
  int p_end = std::min(p_start + rows_per_proc, n) - 1;
  int p_rows = (p_start < n) ? (p_end - p_start + 1) : 0;

  if (p_rows == 0) {
    int zero = 0;
    MPI_Send(&zero, 1, MPI_INT, proc, 0, MPI_COMM_WORLD);
    return;
  }

  int p_halo_first = std::max(0, p_start - 1);
  int p_halo_last = std::min(n - 1, p_end + 1);
  int p_halo_rows = p_halo_last - p_halo_first + 1;

  MPI_Send(&p_halo_rows, 1, MPI_INT, proc, 0, MPI_COMM_WORLD);
  MPI_Send(&p_halo_first, 1, MPI_INT, proc, 1, MPI_COMM_WORLD);
  MPI_Send(&p_start, 1, MPI_INT, proc, 2, MPI_COMM_WORLD);
  MPI_Send(&p_rows, 1, MPI_INT, proc, 3, MPI_COMM_WORLD);

  for (int i = 0; i < p_halo_rows; ++i) {
    const auto &row = global_data[p_halo_first + i];
    MPI_Send(row.data(), m, MPI_INT, proc, 4 + i, MPI_COMM_WORLD);
  }
}

void BoltenkovSGaussianKernelALL::FillLocalHaloForRoot(int local_start_row, int local_end_row, int n, int m,
                                                       const std::vector<std::vector<int>> &global_data,
                                                       std::vector<std::vector<int>> &local_halo) {
  int halo_first = std::max(0, local_start_row - 1);
  int halo_last = std::min(n - 1, local_end_row + 1);
  int halo_rows = halo_last - halo_first + 1;

  if (halo_rows > 0 && m > 0 && halo_rows < 1e6 && m < 1e6) {
    local_halo.resize(halo_rows);
    for (int i = 0; i < halo_rows; ++i) {
      local_halo[i].resize(m);
    }
  } else {
    return;
  }

  for (int i = 0; i < halo_rows; ++i) {
    local_halo[i] = global_data[halo_first + i];
  }
}

void BoltenkovSGaussianKernelALL::ReceiveRowsOnWorker(int m, int &local_start_row, int &local_rows,
                                                      std::vector<std::vector<int>> &local_halo) {
  int recv_halo_rows = 0;
  MPI_Recv(&recv_halo_rows, 1, MPI_INT, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

  if (recv_halo_rows == 0) {
    local_rows = 0;
    return;
  }

  int recv_halo_first = 0;
  MPI_Recv(&recv_halo_first, 1, MPI_INT, 0, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
  MPI_Recv(&local_start_row, 1, MPI_INT, 0, 2, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
  MPI_Recv(&local_rows, 1, MPI_INT, 0, 3, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

  if (recv_halo_rows > 0 && m > 0 && recv_halo_rows < 1e6 && m < 1e6) {
    local_halo.resize(recv_halo_rows);
    for (int i = 0; i < recv_halo_rows; ++i) {
      local_halo[i].resize(m);
    }
  } else {
    return;
  }
  for (int i = 0; i < recv_halo_rows; ++i) {
    MPI_Recv(local_halo[i].data(), m, MPI_INT, 0, 4 + i, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
  }
}

void BoltenkovSGaussianKernelALL::ScatterRows(std::vector<std::vector<int>> &global_data,
                                              std::vector<std::vector<int>> &local_halo, int &local_start_row,
                                              int &local_rows, int m, int rank, int size) {
  int n = static_cast<int>(global_data.size());
  int rows_per_proc = (n + size - 1) / size;

  local_start_row = rank * rows_per_proc;
  int local_end_row = std::min(local_start_row + rows_per_proc, n) - 1;
  local_rows = (local_start_row < n) ? (local_end_row - local_start_row + 1) : 0;

  if (rank == 0) {
    for (int proc = 1; proc < size; ++proc) {
      SendRowsToOneProcess(proc, rows_per_proc, n, m, global_data);
    }
    if (local_rows > 0) {
      FillLocalHaloForRoot(local_start_row, local_end_row, n, m, global_data, local_halo);
    }
  } else {
    ReceiveRowsOnWorker(m, local_start_row, local_rows, local_halo);
  }
}

void BoltenkovSGaussianKernelALL::GatherResultsRoot(std::vector<std::vector<int>> &local_res, int local_start_row,
                                                    int local_rows, int m, int size) {
  for (int i = 0; i < local_rows; ++i) {
    GetOutput()[local_start_row + i] = local_res[i];
  }

  for (int proc = 1; proc < size; ++proc) {
    int p_rows = 0;
    MPI_Recv(&p_rows, 1, MPI_INT, proc, 10, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    if (p_rows == 0) {
      continue;
    }

    int p_start = 0;
    MPI_Recv(&p_start, 1, MPI_INT, proc, 11, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

    std::vector<std::vector<int>> p_res = std::vector<std::vector<int>>();
    if (p_rows > 0 && m > 0 && p_rows < 1e6 && m < 1e6) {
      p_res.resize(p_rows);
      for (int i = 0; i < p_rows; ++i) {
        p_res[i].resize(m);
      }
    } else {
      return;
    }
    for (int i = 0; i < p_rows; ++i) {
      MPI_Recv(p_res[i].data(), m, MPI_INT, proc, 12 + i, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    }
    for (int i = 0; i < p_rows; ++i) {
      GetOutput()[p_start + i] = p_res[i];
    }
  }
}

void BoltenkovSGaussianKernelALL::GatherResultsOthers(std::vector<std::vector<int>> &local_res, int local_start_row,
                                                      int local_rows, int m) {
  if (local_rows > 0) {
    MPI_Send(&local_rows, 1, MPI_INT, 0, 10, MPI_COMM_WORLD);
    MPI_Send(&local_start_row, 1, MPI_INT, 0, 11, MPI_COMM_WORLD);
    for (int i = 0; i < local_rows; ++i) {
      MPI_Send(local_res[i].data(), m, MPI_INT, 0, 12 + i, MPI_COMM_WORLD);
    }
  } else {
    int zero = 0;
    MPI_Send(&zero, 1, MPI_INT, 0, 10, MPI_COMM_WORLD);
  }
}

void BoltenkovSGaussianKernelALL::GatherResults(std::vector<std::vector<int>> &local_res, int local_start_row,
                                                int local_rows, int m, int rank, int size) {
  if (rank == 0) {
    GatherResultsRoot(local_res, local_start_row, local_rows, m, size);
  } else {
    GatherResultsOthers(local_res, local_start_row, local_rows, m);
  }
}

std::vector<std::vector<int>> BoltenkovSGaussianKernelALL::CreateValidMatrix(int rows, int cols) {
  std::vector<std::vector<int>> result;
  if (rows > 0 && cols > 0 && rows < 1e6 && cols < 1e6) {
    result.resize(rows);
    for (int i = 0; i < rows; ++i) {
      result[i].resize(cols, 0);
    }
  }
  return result;
}

void BoltenkovSGaussianKernelALL::FillTmpFromHalo(std::vector<std::vector<int>> &tmp,
                                                  const std::vector<std::vector<int>> &local_halo, int local_start_row,
                                                  int local_rows, int m) {
  int halo_first = std::max(0, local_start_row - 1);
  int halo_rows = static_cast<int>(local_halo.size());

  for (int i = 0; i < local_rows + 2; ++i) {
    int global_row = local_start_row - 1 + i;
    if (global_row >= halo_first && global_row < halo_first + halo_rows) {
      const auto &src = local_halo[global_row - halo_first];
      auto dst = tmp[i].begin() + 1;
      for (int j = 0; j < m; ++j) {
        dst[j] = src[j];
      }
    }
  }
}

std::vector<std::vector<int>> BoltenkovSGaussianKernelALL::ApplyGaussianFilter(
    const std::vector<std::vector<int>> &local_halo, int local_start_row, int local_rows, int m) {
  auto tmp = CreateValidMatrix(local_rows + 2, m + 2);
  if (tmp.empty()) {
    return {};
  }

  FillTmpFromHalo(tmp, local_halo, local_start_row, local_rows, m);

  auto local_res = CreateValidMatrix(local_rows, m);
  if (local_res.empty()) {
    return {};
  }

  const auto &kernel = kernel_;
  int shift = shift_;

#pragma omp parallel for num_threads(ppc::util::GetNumThreads()) default(none) \
    shared(tmp, local_res, local_rows, m, kernel, shift)
  for (int i = 0; i < local_rows; ++i) {
    for (int j = 1; j <= m; ++j) {
      int val = (tmp[i][j - 1] * kernel[0][0]) + (tmp[i][j] * kernel[0][1]) + (tmp[i][j + 1] * kernel[0][2]) +
                (tmp[i + 1][j - 1] * kernel[1][0]) + (tmp[i + 1][j] * kernel[1][1]) +
                (tmp[i + 1][j + 1] * kernel[1][2]) + (tmp[i + 2][j - 1] * kernel[2][0]) +
                (tmp[i + 2][j] * kernel[2][1]) + (tmp[i + 2][j + 1] * kernel[2][2]);
      local_res[i][j - 1] = val >> shift;
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
  if (n <= 0 || m <= 0 || n >= 1e6 || m >= 1e6) {
    return false;
  }
  BcastSizes(n, m, rank);

  std::vector<std::vector<int>> global_data = std::vector<std::vector<int>>();
  if (rank == 0) {
    global_data = std::get<2>(GetInput());
  }

  std::vector<std::vector<int>> local_halo;
  int local_start_row = 0;
  int local_rows = 0;
  ScatterRows(global_data, local_halo, local_start_row, local_rows, m, rank, size);

  std::vector<std::vector<int>> local_res;
  if (local_rows > 0) {
    local_res = ApplyGaussianFilter(local_halo, local_start_row, local_rows, m);
  }

  GatherResults(local_res, local_start_row, local_rows, m, rank, size);

  if (m > 0) {
    for (int i = 0; i < n; ++i) {
      MPI_Bcast(GetOutput()[i].data(), m, MPI_INT, 0, MPI_COMM_WORLD);
    }
  }

  return true;
}

bool BoltenkovSGaussianKernelALL::PostProcessingImpl() {
  return true;
}

}  // namespace boltenkov_s_gaussian_kernel
