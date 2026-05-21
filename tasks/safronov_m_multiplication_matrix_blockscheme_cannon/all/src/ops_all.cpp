#include "safronov_m_multiplication_matrix_blockscheme_cannon/all/include/ops_all.hpp"

#include <mpi.h>

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <utility>
#include <vector>

#include "oneapi/tbb/blocked_range2d.h"
#include "oneapi/tbb/parallel_for.h"
#include "safronov_m_multiplication_matrix_blockscheme_cannon/common/include/common.hpp"

namespace safronov_m_multiplication_matrix_blocksscheme_cannon {

SafronovMMultiplicationMatrixBlockSchemeCannonALL::SafronovMMultiplicationMatrixBlockSchemeCannonALL(const InType &in) {
  SetTypeOfTask(GetStaticTypeOfTask());
  GetInput() = in;
}

bool SafronovMMultiplicationMatrixBlockSchemeCannonALL::ValidationImpl() {
  const auto &in = GetInput();
  int n = std::get<0>(in);
  const auto &a = std::get<1>(in);
  const auto &b = std::get<2>(in);

  return (n > 0) && (!a.empty() && !b.empty()) && (a.size() == a[0].size()) && (b.size() == b[0].size()) &&
         (a.size() == b.size());
}

bool SafronovMMultiplicationMatrixBlockSchemeCannonALL::PreProcessingImpl() {
  GetOutput().clear();
  return true;
}

int SafronovMMultiplicationMatrixBlockSchemeCannonALL::CalcPaddedSize(int n, int q) {
  return ((n + q - 1) / q) * q;
}

void SafronovMMultiplicationMatrixBlockSchemeCannonALL::PadMatrix(const std::vector<std::vector<double>> &src,
                                                                  std::vector<std::vector<double>> &dst, int padded_n) {
  auto sz = static_cast<size_t>(padded_n);

  dst = std::vector<std::vector<double>>(sz, std::vector<double>(sz, 0.0));

  size_t n = src.size();
  for (size_t i = 0; i < n; ++i) {
    std::copy(src[i].begin(), src[i].end(), dst[i].begin());
  }
}

void SafronovMMultiplicationMatrixBlockSchemeCannonALL::ParallelMultiplyBlocks(const std::vector<double> &a,
                                                                               const std::vector<double> &b,
                                                                               std::vector<double> &c, int block_size) {
  tbb::parallel_for(tbb::blocked_range2d<int>(0, block_size, 0, block_size), [&](const tbb::blocked_range2d<int> &r) {
    for (int i = r.rows().begin(); i < r.rows().end(); ++i) {
      for (int k = 0; k < block_size; ++k) {
        double temp = a[(i * block_size) + k];

        for (int j = r.cols().begin(); j < r.cols().end(); ++j) {
          c[(i * block_size) + j] += temp * b[(k * block_size) + j];
        }
      }
    }
  });
}

void SafronovMMultiplicationMatrixBlockSchemeCannonALL::DistributeData(
    MPI_Comm comm, int worker_rank, int worker_size, int q, int block_size,
    const std::vector<std::vector<double>> &matrix_a_full, const std::vector<std::vector<double>> &matrix_b_full,
    std::vector<double> &local_a, std::vector<double> &local_b) {
  size_t b_sz = static_cast<size_t>(block_size) * block_size;
  if (worker_rank == 0) {
    for (int proc = 0; proc < worker_size; ++proc) {
      int row = proc / q;
      int col = proc % q;

      std::vector<double> send_a(static_cast<size_t>(block_size) * block_size);
      std::vector<double> send_b(static_cast<size_t>(block_size) * block_size);

      for (int i = 0; i < block_size; ++i) {
        for (int j = 0; j < block_size; ++j) {
          int a_row = (row * block_size) + i;
          int a_col = (((col + row) % q) * block_size) + j;
          int b_row = (((row + col) % q) * block_size) + i;
          int b_col = (col * block_size) + j;

          send_a[(i * block_size) + j] = matrix_a_full[a_row][a_col];
          send_b[(i * block_size) + j] = matrix_b_full[b_row][b_col];
        }
      }

      if (proc == 0) {
        local_a = std::move(send_a);
        local_b = std::move(send_b);
      } else {
        MPI_Send(send_a.data(), block_size * block_size, MPI_DOUBLE, proc, 0, comm);
        MPI_Send(send_b.data(), block_size * block_size, MPI_DOUBLE, proc, 1, comm);
      }
    }
  } else {
    local_a.resize(b_sz);
    local_b.resize(b_sz);
    MPI_Recv(local_a.data(), block_size * block_size, MPI_DOUBLE, 0, 0, comm, MPI_STATUS_IGNORE);
    MPI_Recv(local_b.data(), block_size * block_size, MPI_DOUBLE, 0, 1, comm, MPI_STATUS_IGNORE);
  }
}

void SafronovMMultiplicationMatrixBlockSchemeCannonALL::CannonAlgorithm(MPI_Comm comm, int worker_rank, int q,
                                                                        int block_size, std::vector<double> &local_a,
                                                                        std::vector<double> &local_b,
                                                                        std::vector<double> &local_c) {
  int row = worker_rank / q;
  int col = worker_rank % q;

  int left = (row * q) + ((col - 1 + q) % q);
  int right = (row * q) + ((col + 1) % q);
  int up = (((row - 1 + q) % q) * q) + col;
  int down = (((row + 1) % q) * q) + col;

  for (int step = 0; step < q; ++step) {
    ParallelMultiplyBlocks(local_a, local_b, local_c, block_size);

    if (step < q - 1) {
      std::vector<double> next_a(static_cast<size_t>(block_size) * block_size);
      std::vector<double> next_b(static_cast<size_t>(block_size) * block_size);

      MPI_Sendrecv(local_a.data(), block_size * block_size, MPI_DOUBLE, left, 10, next_a.data(),
                   block_size * block_size, MPI_DOUBLE, right, 10, comm, MPI_STATUS_IGNORE);

      MPI_Sendrecv(local_b.data(), block_size * block_size, MPI_DOUBLE, up, 11, next_b.data(), block_size * block_size,
                   MPI_DOUBLE, down, 11, comm, MPI_STATUS_IGNORE);

      local_a = std::move(next_a);
      local_b = std::move(next_b);
    }
  }
}

void SafronovMMultiplicationMatrixBlockSchemeCannonALL::FillResultFromBuffer(std::vector<double> &flat_result,
                                                                             const std::vector<double> &buffer, int row,
                                                                             int col, int block_size, int padded_n) {
  for (int i = 0; i < block_size; ++i) {
    for (int j = 0; j < block_size; ++j) {
      int global_row = (row * block_size) + i;
      int global_col = (col * block_size) + j;
      flat_result[(global_row * padded_n) + global_col] = buffer[(i * block_size) + j];
    }
  }
}

void SafronovMMultiplicationMatrixBlockSchemeCannonALL::CollectResult(MPI_Comm comm, int worker_rank, int worker_size,
                                                                      int q, int block_size,
                                                                      std::vector<double> &flat_result,
                                                                      const std::vector<double> &local_c) {
  int padded_n = q * block_size;

  if (worker_rank == 0) {
    FillResultFromBuffer(flat_result, local_c, 0, 0, block_size, padded_n);

    std::vector<double> recv_buf(static_cast<std::size_t>(block_size) * block_size);
    for (int proc = 1; proc < worker_size; ++proc) {
      MPI_Recv(recv_buf.data(), block_size * block_size, MPI_DOUBLE, proc, 20, comm, MPI_STATUS_IGNORE);
      FillResultFromBuffer(flat_result, recv_buf, proc / q, proc % q, block_size, padded_n);
    }
  } else {
    MPI_Send(local_c.data(), block_size * block_size, MPI_DOUBLE, 0, 20, comm);
  }
}

bool SafronovMMultiplicationMatrixBlockSchemeCannonALL::RunImpl() {
  int rank = 0;
  int size = 1;

  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &size);

  int q = static_cast<int>(std::sqrt(size));
  int active = q * q;

  int n = 0;
  if (rank == 0) {
    n = static_cast<int>(std::get<1>(GetInput()).size());
  }

  MPI_Bcast(&n, 1, MPI_INT, 0, MPI_COMM_WORLD);

  int padded_n = CalcPaddedSize(n, std::max(1, q));
  int block_size = padded_n / std::max(1, q);

  std::vector<std::vector<double>> padded_a;
  std::vector<std::vector<double>> padded_b;

  if (rank == 0) {
    PadMatrix(std::get<1>(GetInput()), padded_a, padded_n);
    PadMatrix(std::get<2>(GetInput()), padded_b, padded_n);
  }

  std::vector<double> flat_result(static_cast<size_t>(padded_n) * padded_n, 0.0);

  MPI_Comm comm = MPI_COMM_NULL;
  int color = (rank < active) ? 0 : MPI_UNDEFINED;

  MPI_Comm_split(MPI_COMM_WORLD, color, rank, &comm);

  if (rank < active) {
    int wrank = 0;
    int wsize = 0;

    MPI_Comm_rank(comm, &wrank);
    MPI_Comm_size(comm, &wsize);

    std::vector<double> local_a(static_cast<size_t>(block_size) * block_size);
    std::vector<double> local_b(static_cast<size_t>(block_size) * block_size);
    std::vector<double> local_c(static_cast<size_t>(block_size) * block_size, 0.0);

    DistributeData(comm, wrank, wsize, q, block_size, padded_a, padded_b, local_a, local_b);

    CannonAlgorithm(comm, wrank, q, block_size, local_a, local_b, local_c);

    CollectResult(comm, wrank, wsize, q, block_size, flat_result, local_c);

    MPI_Comm_free(&comm);
  }

  MPI_Bcast(flat_result.data(), padded_n * padded_n, MPI_DOUBLE, 0, MPI_COMM_WORLD);

  std::vector<std::vector<double>> result(static_cast<size_t>(n), std::vector<double>(static_cast<size_t>(n)));

  auto un = static_cast<size_t>(n);
  auto upadded_n = static_cast<size_t>(padded_n);

  for (size_t i = 0; i < un; ++i) {
    for (size_t j = 0; j < un; ++j) {
      result[i][j] = flat_result[(i * upadded_n) + j];
    }
  }

  GetOutput() = std::move(result);
  return true;
}

bool SafronovMMultiplicationMatrixBlockSchemeCannonALL::PostProcessingImpl() {
  return true;
}

}  // namespace safronov_m_multiplication_matrix_blocksscheme_cannon
