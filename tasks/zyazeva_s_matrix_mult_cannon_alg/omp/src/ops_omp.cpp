#include "zyazeva_s_matrix_mult_cannon_alg/omp/include/ops_omp.hpp"

#include <omp.h>

#include <cmath>
#include <cstddef>
#include <utility>
#include <vector>

#include "zyazeva_s_matrix_mult_cannon_alg/common/include/common.hpp"

namespace zyazeva_s_matrix_mult_cannon_alg {

bool ZyazevaSMatrixMultCannonAlgOMP::IsPerfectSquare(int x) {
  int root = static_cast<int>(std::sqrt(x));
  return root * root == x;
}

void ZyazevaSMatrixMultCannonAlgOMP::MultiplyBlocks(const std::vector<double> &a, const std::vector<double> &b,
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

ZyazevaSMatrixMultCannonAlgOMP::ZyazevaSMatrixMultCannonAlgOMP(const InType &in) {
  SetTypeOfTask(GetStaticTypeOfTask());
  GetInput() = in;
  GetOutput() = {};
}

bool ZyazevaSMatrixMultCannonAlgOMP::ValidationImpl() {
  const size_t sz = std::get<0>(GetInput());
  const auto &m1 = std::get<1>(GetInput());
  const auto &m2 = std::get<2>(GetInput());

  return sz > 0 && m1.size() == sz * sz && m2.size() == sz * sz;
}

bool ZyazevaSMatrixMultCannonAlgOMP::PreProcessingImpl() {
  GetOutput() = {};
  return true;
}

void ZyazevaSMatrixMultCannonAlgOMP::RegularMultiplication(const std::vector<double> &m1, const std::vector<double> &m2,
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

void ZyazevaSMatrixMultCannonAlgOMP::InitializeBlocks(const std::vector<double> &m1, const std::vector<double> &m2,
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

void ZyazevaSMatrixMultCannonAlgOMP::AlignBlocks(const std::vector<std::vector<double>> &blocks_a,
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

void ZyazevaSMatrixMultCannonAlgOMP::CannonStep(std::vector<std::vector<double>> &aligned_a,
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

void ZyazevaSMatrixMultCannonAlgOMP::AssembleResult(const std::vector<std::vector<double>> &blocks_c,
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

bool ZyazevaSMatrixMultCannonAlgOMP::RunImpl() {
  const auto sz = static_cast<int>(std::get<0>(GetInput()));
  const auto &m1 = std::get<1>(GetInput());
  const auto &m2 = std::get<2>(GetInput());

  std::vector<double> res_m(static_cast<size_t>(sz) * static_cast<size_t>(sz), 0.0);

  int num_threads = 1;
#pragma omp parallel default(none) shared(num_threads)
  {
#pragma omp single
    num_threads = omp_get_num_threads();
  }

  const bool can_use_cannon =
      IsPerfectSquare(num_threads) && sz >= num_threads && (sz % static_cast<int>(std::sqrt(num_threads)) == 0);

  if (!can_use_cannon) {
    RegularMultiplication(m1, m2, res_m, sz);
    GetOutput() = res_m;
    return true;
  }

  const int grid_size = static_cast<int>(std::sqrt(num_threads));
  const int block_size = sz / grid_size;

  const auto grid_size_t = static_cast<size_t>(grid_size);
  const auto block_size_t = static_cast<size_t>(block_size);
  const auto sz_t = static_cast<size_t>(sz);

  std::vector<std::vector<double>> blocks_a(grid_size_t * grid_size_t);
  std::vector<std::vector<double>> blocks_b(grid_size_t * grid_size_t);
  std::vector<std::vector<double>> blocks_c(grid_size_t * grid_size_t,
                                            std::vector<double>(block_size_t * block_size_t, 0.0));

  InitializeBlocks(m1, m2, blocks_a, blocks_b, grid_size, block_size, grid_size_t, block_size_t, sz_t);

  std::vector<std::vector<double>> aligned_a(grid_size_t * grid_size_t);
  std::vector<std::vector<double>> aligned_b(grid_size_t * grid_size_t);
  AlignBlocks(blocks_a, blocks_b, aligned_a, aligned_b, grid_size, grid_size_t);

  for (int step = 0; step < grid_size; ++step) {
    CannonStep(aligned_a, aligned_b, blocks_c, grid_size, block_size, grid_size_t, step);
  }

  AssembleResult(blocks_c, res_m, grid_size, block_size, sz_t, grid_size_t, block_size_t);

  GetOutput() = res_m;
  return true;
}

bool ZyazevaSMatrixMultCannonAlgOMP::PostProcessingImpl() {
  return true;
}

}  // namespace zyazeva_s_matrix_mult_cannon_alg
