#include "chyokotov_a_dense_matrix_mul_foxs_algorithm/all/include/ops_all.hpp"

#include <mpi.h>
#include <tbb/blocked_range2d.h>
#include <tbb/parallel_for.h>

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <utility>
#include <vector>

#include "chyokotov_a_dense_matrix_mul_foxs_algorithm/common/include/common.hpp"

namespace chyokotov_a_dense_matrix_mul_foxs_algorithm {

ChyokotovADenseMatMulFoxAlgorithmALL::ChyokotovADenseMatMulFoxAlgorithmALL(const InType &in) {
  SetTypeOfTask(GetStaticTypeOfTask());
  GetInput() = in;
  GetOutput().clear();
}

bool ChyokotovADenseMatMulFoxAlgorithmALL::ValidationImpl() {
  return (GetInput().first.size() == GetInput().second.size());
}

bool ChyokotovADenseMatMulFoxAlgorithmALL::PreProcessingImpl() {
  GetOutput().clear();
  GetOutput().resize(GetInput().first.size(), 0.0);
  return true;
}

int ChyokotovADenseMatMulFoxAlgorithmALL::CalcPaddedSize(int n, int q) {
  if (q <= 0) {
    return n;
  }
  return ((n + q - 1) / q) * q;
}

void ChyokotovADenseMatMulFoxAlgorithmALL::PadMatrix(const std::vector<double> &src, std::vector<double> &dst,
                                                     int original_n, int padded_n) {
  dst.assign(static_cast<size_t>(padded_n) * padded_n, 0.0);

  for (int i = 0; i < original_n; ++i) {
    for (int j = 0; j < original_n; ++j) {
      dst[(i * padded_n) + j] = src[(i * original_n) + j];
    }
  }
}

void ChyokotovADenseMatMulFoxAlgorithmALL::Multiply(const std::vector<double> &a_block,
                                                    const std::vector<double> &b_block, std::vector<double> &c_block,
                                                    int block_size) {
  tbb::parallel_for(tbb::blocked_range2d<int>(0, block_size, 0, block_size),
                    [&](const tbb::blocked_range2d<int> &range) {
    for (int i = range.rows().begin(); i < range.rows().end(); ++i) {
      for (int k = 0; k < block_size; ++k) {
        double temp = a_block[(i * block_size) + k];
        for (int j = range.cols().begin(); j < range.cols().end(); ++j) {
          c_block[(i * block_size) + j] += temp * b_block[(k * block_size) + j];
        }
      }
    }
  });
}

void ChyokotovADenseMatMulFoxAlgorithmALL::DistributeData(MPI_Comm comm, int worker_rank, int worker_size, int q,
                                                          int block_size, const std::vector<double> &matrix_a_full,
                                                          const std::vector<double> &matrix_b_full,
                                                          std::vector<double> &local_a, std::vector<double> &local_b) {
  size_t block_sz = static_cast<size_t>(block_size) * block_size;

  if (worker_rank == 0) {
    for (int proc = 0; proc < worker_size; ++proc) {
      int row = proc / q;
      int col = proc % q;

      std::vector<double> send_a(block_sz);
      std::vector<double> send_b(block_sz);

      for (int i = 0; i < block_size; ++i) {
        for (int j = 0; j < block_size; ++j) {
          int a_row = (row * block_size) + i;
          int a_col = (((col + row) % q) * block_size) + j;
          int b_row = (((row + col) % q) * block_size) + i;
          int b_col = (col * block_size) + j;

          send_a[(i * block_size) + j] = matrix_a_full[(a_row * block_size * q) + a_col];
          send_b[(i * block_size) + j] = matrix_b_full[(b_row * block_size * q) + b_col];
        }
      }

      if (proc == 0) {
        local_a = std::move(send_a);
        local_b = std::move(send_b);
      } else {
        MPI_Send(send_a.data(), static_cast<int>(block_sz), MPI_DOUBLE, proc, 0, comm);
        MPI_Send(send_b.data(), static_cast<int>(block_sz), MPI_DOUBLE, proc, 1, comm);
      }
    }
  } else {
    local_a.resize(block_sz);
    local_b.resize(block_sz);
    MPI_Recv(local_a.data(), static_cast<int>(block_sz), MPI_DOUBLE, 0, 0, comm, MPI_STATUS_IGNORE);
    MPI_Recv(local_b.data(), static_cast<int>(block_sz), MPI_DOUBLE, 0, 1, comm, MPI_STATUS_IGNORE);
  }
}

void ChyokotovADenseMatMulFoxAlgorithmALL::FoxAlgorithm(MPI_Comm comm, int worker_rank, int q, int block_size,
                                                        std::vector<double> &local_a, std::vector<double> &local_b,
                                                        std::vector<double> &local_c) {
  int row = worker_rank / q;
  int col = worker_rank % q;

  int left = (row * q) + ((col - 1 + q) % q);
  int right = (row * q) + ((col + 1) % q);
  int up = (((row - 1 + q) % q) * q) + col;
  int down = (((row + 1) % q) * q) + col;

  size_t block_sz = static_cast<size_t>(block_size) * block_size;
  std::vector<double> next_a(block_sz);
  std::vector<double> next_b(block_sz);

  for (int step = 0; step < q; ++step) {
    Multiply(local_a, local_b, local_c, block_size);

    if (step < q - 1) {
      MPI_Sendrecv(local_a.data(), static_cast<int>(block_sz), MPI_DOUBLE, left, 10, next_a.data(),
                   static_cast<int>(block_sz), MPI_DOUBLE, right, 10, comm, MPI_STATUS_IGNORE);

      MPI_Sendrecv(local_b.data(), static_cast<int>(block_sz), MPI_DOUBLE, up, 11, next_b.data(),
                   static_cast<int>(block_sz), MPI_DOUBLE, down, 11, comm, MPI_STATUS_IGNORE);

      local_a.swap(next_a);
      local_b.swap(next_b);
    }
  }
}

void ChyokotovADenseMatMulFoxAlgorithmALL::CollectResult(MPI_Comm comm, int worker_rank, int worker_size, int q,
                                                         int block_size, std::vector<double> &flat_result,
                                                         const std::vector<double> &local_c) {
  int padded_n = q * block_size;

  auto fillres = [&](const std::vector<double> &buffer, int row, int col) {
    for (int i = 0; i < block_size; ++i) {
      for (int j = 0; j < block_size; ++j) {
        int global_row = (row * block_size) + i;
        int global_col = (col * block_size) + j;
        flat_result[(global_row * padded_n) + global_col] = buffer[(i * block_size) + j];
      }
    }
  };

  if (worker_rank == 0) {
    fillres(local_c, 0, 0);

    std::vector<double> recv_buf(static_cast<size_t>(block_size) * block_size);
    for (int proc = 1; proc < worker_size; ++proc) {
      MPI_Recv(recv_buf.data(), static_cast<int>(recv_buf.size()), MPI_DOUBLE, proc, 20, comm, MPI_STATUS_IGNORE);
      fillres(recv_buf, proc / q, proc % q);
    }
  } else {
    MPI_Send(local_c.data(), static_cast<int>(local_c.size()), MPI_DOUBLE, 0, 20, comm);
  }
}

bool ChyokotovADenseMatMulFoxAlgorithmALL::RunImpl() {
  int rank = 0;
  int size = 1;

  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &size);

  int q = static_cast<int>(std::sqrt(static_cast<double>(size)));
  int active = q * q;

  int n = 0;
  if (rank == 0) {
    n = static_cast<int>(std::sqrt(static_cast<double>(GetInput().first.size())));
  }

  MPI_Bcast(&n, 1, MPI_INT, 0, MPI_COMM_WORLD);

  if (n == 0) {
    return true;
  }

  int padded_n = CalcPaddedSize(n, std::max(1, q));
  int block_size = padded_n / std::max(1, q);

  std::vector<double> padded_a;
  std::vector<double> padded_b;

  if (rank == 0) {
    PadMatrix(GetInput().first, padded_a, n, padded_n);
    PadMatrix(GetInput().second, padded_b, n, padded_n);
  }

  MPI_Comm comm = MPI_COMM_NULL;
  int color = (rank < active) ? 0 : MPI_UNDEFINED;
  MPI_Comm_split(MPI_COMM_WORLD, color, rank, &comm);

  std::vector<double> flat_result(static_cast<size_t>(padded_n) * padded_n, 0.0);

  if (rank < active) {
    int wrank = 0;
    int wsize = 0;
    MPI_Comm_rank(comm, &wrank);
    MPI_Comm_size(comm, &wsize);

    size_t block_sz = static_cast<size_t>(block_size) * block_size;
    std::vector<double> local_a(block_sz);
    std::vector<double> local_b(block_sz);
    std::vector<double> local_c(block_sz, 0.0);

    DistributeData(comm, wrank, wsize, q, block_size, padded_a, padded_b, local_a, local_b);
    FoxAlgorithm(comm, wrank, q, block_size, local_a, local_b, local_c);
    CollectResult(comm, wrank, wsize, q, block_size, flat_result, local_c);

    MPI_Comm_free(&comm);
  }

  MPI_Bcast(flat_result.data(), padded_n * padded_n, MPI_DOUBLE, 0, MPI_COMM_WORLD);

  GetOutput().resize(static_cast<size_t>(n) * n);
  for (int i = 0; i < n; ++i) {
    for (int j = 0; j < n; ++j) {
      GetOutput()[(i * n) + j] = flat_result[(i * padded_n) + j];
    }
  }

  return true;
}

bool ChyokotovADenseMatMulFoxAlgorithmALL::PostProcessingImpl() {
  return true;
}

}  // namespace chyokotov_a_dense_matrix_mul_foxs_algorithm
