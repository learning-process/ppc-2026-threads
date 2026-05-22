#include "zyazeva_s_matrix_mult_cannon_alg/all/include/ops_all.hpp"

#include <mpi.h>
#include <omp.h>

#include <cmath>
#include <cstddef>
#include <utility>
#include <vector>

#include "zyazeva_s_matrix_mult_cannon_alg/common/include/common.hpp"

namespace zyazeva_s_matrix_mult_cannon_alg {

bool ZyazevaSMatrixMultCannonAlgALL::IsPerfectSquare(int x) {
  int root = static_cast<int>(std::sqrt(x));
  return root * root == x;
}

void ZyazevaSMatrixMultCannonAlgALL::MultiplyBlocks(const std::vector<double> &a, const std::vector<double> &b,
                                                    std::vector<double> &c, int block_size) {
  for (int i = 0; i < block_size; ++i) {
    for (int k = 0; k < block_size; ++k) {
      const size_t i_idx = static_cast<size_t>(i) * static_cast<size_t>(block_size);
      const size_t k_idx = static_cast<size_t>(k) * static_cast<size_t>(block_size);
      double a_ik = a[i_idx + static_cast<size_t>(k)];
      for (int j = 0; j < block_size; ++j) {
        c[i_idx + static_cast<size_t>(j)] += a_ik * b[k_idx + static_cast<size_t>(j)];
      }
    }
  }
}

ZyazevaSMatrixMultCannonAlgALL::ZyazevaSMatrixMultCannonAlgALL(const InType &in) {
  SetTypeOfTask(GetStaticTypeOfTask());
  GetInput() = in;
  GetOutput() = {};
  MPI_Comm_rank(MPI_COMM_WORLD, &rank_);
  MPI_Comm_size(MPI_COMM_WORLD, &mpi_size_);
}

bool ZyazevaSMatrixMultCannonAlgALL::ValidationImpl() {
  int valid = 0;
  if (rank_ == 0) {
    const size_t sz = std::get<0>(GetInput());
    const auto &m1 = std::get<1>(GetInput());
    const auto &m2 = std::get<2>(GetInput());
    valid = (sz > 0 && m1.size() == sz * sz && m2.size() == sz * sz) ? 1 : 0;
  }
  MPI_Bcast(&valid, 1, MPI_INT, 0, MPI_COMM_WORLD);
  return valid != 0;
}

bool ZyazevaSMatrixMultCannonAlgALL::PreProcessingImpl() {
  GetOutput() = {};
  return true;
}

void ZyazevaSMatrixMultCannonAlgALL::RegularMultiplication(const std::vector<double> &m1, const std::vector<double> &m2,
                                                           std::vector<double> &res, int sz) {
#pragma omp parallel for default(none) shared(m1, m2, res, sz)
  for (int i = 0; i < sz; ++i) {
    const size_t i_offset = static_cast<size_t>(i) * static_cast<size_t>(sz);
    for (int j = 0; j < sz; ++j) {
      double sum = 0.0;
      for (int k = 0; k < sz; ++k) {
        const size_t k_offset = static_cast<size_t>(k) * static_cast<size_t>(sz);
        sum += m1[i_offset + static_cast<size_t>(k)] * m2[k_offset + static_cast<size_t>(j)];
      }
      res[i_offset + static_cast<size_t>(j)] = sum;
    }
  }
}

void ZyazevaSMatrixMultCannonAlgALL::InitializeBlocks(const std::vector<double> &m1, const std::vector<double> &m2,
                                                      std::vector<std::vector<double>> &blocks_a,
                                                      std::vector<std::vector<double>> &blocks_b, int grid_size,
                                                      int block_size, size_t grid_size_t, size_t block_size_t,
                                                      size_t sz_t) {
  for (int i = 0; i < grid_size; ++i) {
    for (int j = 0; j < grid_size; ++j) {
      const size_t block_idx = (static_cast<size_t>(i) * grid_size_t) + static_cast<size_t>(j);
      blocks_a[block_idx].resize(block_size_t * block_size_t);
      blocks_b[block_idx].resize(block_size_t * block_size_t);

      for (int bi = 0; bi < block_size; ++bi) {
        for (int bj = 0; bj < block_size; ++bj) {
          const size_t global_i = (static_cast<size_t>(i) * block_size_t) + static_cast<size_t>(bi);
          const size_t global_j = (static_cast<size_t>(j) * block_size_t) + static_cast<size_t>(bj);
          const size_t local_idx = (static_cast<size_t>(bi) * block_size_t) + static_cast<size_t>(bj);

          blocks_a[block_idx][local_idx] = m1[(global_i * sz_t) + global_j];
          blocks_b[block_idx][local_idx] = m2[(global_i * sz_t) + global_j];
        }
      }
    }
  }
}

void ZyazevaSMatrixMultCannonAlgALL::AlignBlocks(const std::vector<std::vector<double>> &blocks_a,
                                                 const std::vector<std::vector<double>> &blocks_b,
                                                 std::vector<std::vector<double>> &aligned_a,
                                                 std::vector<std::vector<double>> &aligned_b, int grid_size,
                                                 size_t grid_size_t) {
#pragma omp parallel for default(none) shared(blocks_a, blocks_b, aligned_a, aligned_b, grid_size, grid_size_t) \
    collapse(2)
  for (int i = 0; i < grid_size; ++i) {
    for (int j = 0; j < grid_size; ++j) {
      const size_t block_idx = (static_cast<size_t>(i) * grid_size_t) + static_cast<size_t>(j);

      const size_t a_src_idx = (static_cast<size_t>(i) * grid_size_t) + static_cast<size_t>((j + i) % grid_size);
      aligned_a[block_idx] = blocks_a[a_src_idx];

      const size_t b_src_idx = (static_cast<size_t>((i + j) % grid_size) * grid_size_t) + static_cast<size_t>(j);
      aligned_b[block_idx] = blocks_b[b_src_idx];
    }
  }
}

void ZyazevaSMatrixMultCannonAlgALL::CannonStep(std::vector<std::vector<double>> &aligned_a,
                                                std::vector<std::vector<double>> &aligned_b,
                                                std::vector<std::vector<double>> &blocks_c, int grid_size,
                                                int block_size, size_t grid_size_t, int step) {
#pragma omp parallel for default(none) shared(aligned_a, aligned_b, blocks_c, grid_size, block_size, grid_size_t) \
    collapse(2)
  for (int i = 0; i < grid_size; ++i) {
    for (int j = 0; j < grid_size; ++j) {
      const size_t block_idx = (static_cast<size_t>(i) * grid_size_t) + static_cast<size_t>(j);
      MultiplyBlocks(aligned_a[block_idx], aligned_b[block_idx], blocks_c[block_idx], block_size);
    }
  }

  if (step < grid_size - 1) {
    std::vector<std::vector<double>> new_aligned_a(grid_size_t * grid_size_t);
    std::vector<std::vector<double>> new_aligned_b(grid_size_t * grid_size_t);

#pragma omp parallel for default(none) \
    shared(aligned_a, aligned_b, new_aligned_a, new_aligned_b, grid_size, grid_size_t) collapse(2)
    for (int i = 0; i < grid_size; ++i) {
      for (int j = 0; j < grid_size; ++j) {
        const size_t block_idx = (static_cast<size_t>(i) * grid_size_t) + static_cast<size_t>(j);

        const size_t a_src_idx = (static_cast<size_t>(i) * grid_size_t) + static_cast<size_t>((j + 1) % grid_size);
        new_aligned_a[block_idx] = aligned_a[a_src_idx];

        const size_t b_src_idx = (static_cast<size_t>((i + 1) % grid_size) * grid_size_t) + static_cast<size_t>(j);
        new_aligned_b[block_idx] = aligned_b[b_src_idx];
      }
    }

    aligned_a = std::move(new_aligned_a);
    aligned_b = std::move(new_aligned_b);
  }
}

void ZyazevaSMatrixMultCannonAlgALL::AssembleResult(const std::vector<std::vector<double>> &blocks_c,
                                                    std::vector<double> &res_m, int grid_size, int block_size,
                                                    size_t sz_t, size_t grid_size_t, size_t block_size_t) {
#pragma omp parallel for default(none) shared(blocks_c, res_m, grid_size, block_size, sz_t, grid_size_t, block_size_t) \
    collapse(2)
  for (int i = 0; i < grid_size; ++i) {
    for (int j = 0; j < grid_size; ++j) {
      const size_t block_idx = (static_cast<size_t>(i) * grid_size_t) + static_cast<size_t>(j);
      const auto &block = blocks_c[block_idx];

      for (int bi = 0; bi < block_size; ++bi) {
        for (int bj = 0; bj < block_size; ++bj) {
          const size_t global_i = (static_cast<size_t>(i) * block_size_t) + static_cast<size_t>(bi);
          const size_t global_j = (static_cast<size_t>(j) * block_size_t) + static_cast<size_t>(bj);
          const size_t local_idx = (static_cast<size_t>(bi) * block_size_t) + static_cast<size_t>(bj);

          res_m[(global_i * sz_t) + global_j] = block[local_idx];
        }
      }
    }
  }
}

// Extracted from RunImpl to reduce cognitive complexity
void ZyazevaSMatrixMultCannonAlgALL::DistributeBlocks(const std::vector<double> &m1, const std::vector<double> &m2,
                                                      int grid, int block_size, int block_elems, int sz,
                                                      std::vector<double> &local_a,
                                                      std::vector<double> &local_b) const {
  if (rank_ == 0) {
    for (int proc = 0; proc < mpi_size_; ++proc) {
      const int proc_row = proc / grid;
      const int proc_col = proc % grid;

      std::vector<double> tmp_a(block_elems);
      std::vector<double> tmp_b(block_elems);

      for (int i = 0; i < block_size; ++i) {
        for (int j = 0; j < block_size; ++j) {
          const int gi = (proc_row * block_size) + i;
          const int gj = (proc_col * block_size) + j;
          tmp_a[(i * block_size) + j] = m1[(gi * sz) + gj];
          tmp_b[(i * block_size) + j] = m2[(gi * sz) + gj];
        }
      }

      if (proc == 0) {
        local_a = tmp_a;
        local_b = tmp_b;
      } else {
        MPI_Send(tmp_a.data(), block_elems, MPI_DOUBLE, proc, 0, MPI_COMM_WORLD);
        MPI_Send(tmp_b.data(), block_elems, MPI_DOUBLE, proc, 1, MPI_COMM_WORLD);
      }
    }
  } else {
    MPI_Recv(local_a.data(), block_elems, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    MPI_Recv(local_b.data(), block_elems, MPI_DOUBLE, 0, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
  }
}

void ZyazevaSMatrixMultCannonAlgALL::CollectResult(const std::vector<double> &local_c, std::vector<double> &result,
                                                   int grid, int block_size, int block_elems, int sz) const {
  if (rank_ == 0) {
    for (int proc = 0; proc < mpi_size_; ++proc) {
      const int proc_row_p = proc / grid;
      const int proc_col_p = proc % grid;

      std::vector<double> block(block_elems);
      if (proc == 0) {
        block = local_c;
      } else {
        MPI_Recv(block.data(), block_elems, MPI_DOUBLE, proc, 30, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
      }

      for (int i = 0; i < block_size; ++i) {
        for (int j = 0; j < block_size; ++j) {
          const int gi = (proc_row_p * block_size) + i;
          const int gj = (proc_col_p * block_size) + j;
          result[(gi * sz) + gj] = block[(i * block_size) + j];
        }
      }
    }
  } else {
    MPI_Send(local_c.data(), block_elems, MPI_DOUBLE, 0, 30, MPI_COMM_WORLD);
  }
}

bool ZyazevaSMatrixMultCannonAlgALL::RunImpl() {
  int sz = 0;

  if (rank_ == 0) {
    sz = static_cast<int>(std::get<0>(GetInput()));
  }

  MPI_Bcast(&sz, 1, MPI_INT, 0, MPI_COMM_WORLD);

  const auto sz_t = static_cast<size_t>(sz);
  const size_t mat_size = sz_t * sz_t;

  std::vector<double> m1(mat_size);
  std::vector<double> m2(mat_size);

  if (rank_ == 0) {
    m1 = std::get<1>(GetInput());
    m2 = std::get<2>(GetInput());
  }

  MPI_Bcast(m1.data(), static_cast<int>(mat_size), MPI_DOUBLE, 0, MPI_COMM_WORLD);
  MPI_Bcast(m2.data(), static_cast<int>(mat_size), MPI_DOUBLE, 0, MPI_COMM_WORLD);

  std::vector<double> result(mat_size, 0.0);

  const int grid = static_cast<int>(std::sqrt(mpi_size_));
  const bool use_cannon = (grid * grid == mpi_size_) && (grid > 0) && (sz % grid == 0);

  if (!use_cannon) {
    if (rank_ == 0) {
      RegularMultiplication(m1, m2, result, sz);
    }
    MPI_Bcast(result.data(), static_cast<int>(mat_size), MPI_DOUBLE, 0, MPI_COMM_WORLD);
    GetOutput() = result;
    return true;
  }

  const int block_size = sz / grid;
  const int block_elems = block_size * block_size;

  std::vector<double> local_a(block_elems);
  std::vector<double> local_b(block_elems);
  std::vector<double> local_c(block_elems, 0.0);

  DistributeBlocks(m1, m2, grid, block_size, block_elems, sz, local_a, local_b);

  const int proc_row = rank_ / grid;
  const int proc_col = rank_ % grid;

  for (int sh = 0; sh < proc_row; ++sh) {
    const int send_to = (proc_row * grid) + ((proc_col - 1 + grid) % grid);
    const int recv_from = (proc_row * grid) + ((proc_col + 1) % grid);
    MPI_Sendrecv_replace(local_a.data(), block_elems, MPI_DOUBLE, send_to, 10, recv_from, 10, MPI_COMM_WORLD,
                         MPI_STATUS_IGNORE);
  }

  for (int sh = 0; sh < proc_col; ++sh) {
    const int send_to = (((proc_row - 1 + grid) % grid) * grid) + proc_col;
    const int recv_from = (((proc_row + 1) % grid) * grid) + proc_col;
    MPI_Sendrecv_replace(local_b.data(), block_elems, MPI_DOUBLE, send_to, 11, recv_from, 11, MPI_COMM_WORLD,
                         MPI_STATUS_IGNORE);
  }

  for (int step = 0; step < grid; ++step) {
    MultiplyBlocks(local_a, local_b, local_c, block_size);

    const int left = (proc_row * grid) + ((proc_col - 1 + grid) % grid);
    const int right = (proc_row * grid) + ((proc_col + 1) % grid);
    MPI_Sendrecv_replace(local_a.data(), block_elems, MPI_DOUBLE, left, 20, right, 20, MPI_COMM_WORLD,
                         MPI_STATUS_IGNORE);

    const int up = (((proc_row - 1 + grid) % grid) * grid) + proc_col;
    const int down = (((proc_row + 1) % grid) * grid) + proc_col;
    MPI_Sendrecv_replace(local_b.data(), block_elems, MPI_DOUBLE, up, 21, down, 21, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
  }

  CollectResult(local_c, result, grid, block_size, block_elems, sz);

  MPI_Bcast(result.data(), static_cast<int>(mat_size), MPI_DOUBLE, 0, MPI_COMM_WORLD);
  GetOutput() = result;
  return true;
}

bool ZyazevaSMatrixMultCannonAlgALL::PostProcessingImpl() {
  return true;
}

}  // namespace zyazeva_s_matrix_mult_cannon_alg
