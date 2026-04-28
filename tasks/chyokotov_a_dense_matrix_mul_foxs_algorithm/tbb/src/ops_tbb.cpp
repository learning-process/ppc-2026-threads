#include "chyokotov_a_dense_matrix_mul_foxs_algorithm/tbb/include/ops_tbb.hpp"

#include <tbb/blocked_range2d.h>
#include <tbb/parallel_for.h>

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <vector>

#include "chyokotov_a_dense_matrix_mul_foxs_algorithm/common/include/common.hpp"

namespace chyokotov_a_dense_matrix_mul_foxs_algorithm {

ChyokotovADenseMatMulFoxAlgorithmTBB::ChyokotovADenseMatMulFoxAlgorithmTBB(const InType &in) {
  SetTypeOfTask(GetStaticTypeOfTask());
  GetInput() = in;
  GetOutput().clear();
}

bool ChyokotovADenseMatMulFoxAlgorithmTBB::ValidationImpl() {
  return (GetInput().first.size() == GetInput().second.size());
}

bool ChyokotovADenseMatMulFoxAlgorithmTBB::PreProcessingImpl() {
  GetOutput().clear();
  GetOutput().resize(GetInput().first.size(), 0.0);
  return true;
}

int ChyokotovADenseMatMulFoxAlgorithmTBB::CalculateBlockSize(int n) {
  return static_cast<int>(std::sqrt(static_cast<double>(n)));
}

int ChyokotovADenseMatMulFoxAlgorithmTBB::CountBlock(int n, int size) {
  return (n + size - 1) / size;
}

void ChyokotovADenseMatMulFoxAlgorithmTBB::Matmul(std::vector<double> &a, std::vector<double> &b, int n, int istart,
                                                  int iend, int jstart, int jend, int kstart, int kend) {
  const auto n_sz = static_cast<size_t>(n);
  const auto istart_sz = static_cast<size_t>(istart);
  const auto iend_sz = static_cast<size_t>(iend);

  tbb::parallel_for(istart_sz, iend_sz, [&](size_t i) {
    double *output_row = GetOutput().data() + (i * n_sz);
    const double *a_row = a.data() + (i * n_sz);

    for (int j = jstart; j < jend; ++j) {
      long double sum = 0.0L;
      const double *b_col = b.data() + j;

      for (int k = kstart; k < kend; ++k) {
        const auto k_sz = static_cast<size_t>(k);
        sum += static_cast<long double>(a_row[k_sz]) * static_cast<long double>(b_col[k_sz * n_sz]);
      }
      output_row[j] += static_cast<double>(sum);
    }
  });
}

bool ChyokotovADenseMatMulFoxAlgorithmTBB::RunImpl() {
  std::vector<double> a = GetInput().first;
  std::vector<double> b = GetInput().second;
  int n = static_cast<int>(std::sqrt(static_cast<double>(a.size())));
  if (n == 0) {
    return true;
  }

  int block_size = CalculateBlockSize(n);
  int count_block = CountBlock(n, block_size);

  tbb::parallel_for(tbb::blocked_range2d<int>(0, count_block, 0, count_block), [&](const tbb::blocked_range2d<int> &r) {
    for (int ic = r.rows().begin(); ic < r.rows().end(); ++ic) {
      for (int jc = r.cols().begin(); jc < r.cols().end(); ++jc) {
        for (int kc = 0; kc < count_block; ++kc) {
          int istart = ic * block_size;
          int jstart = jc * block_size;
          int kstart = kc * block_size;

          int iend = std::min(istart + block_size, n);
          int jend = std::min(jstart + block_size, n);
          int kend = std::min(kstart + block_size, n);

          Matmul(a, b, n, istart, iend, jstart, jend, kstart, kend);
        }
      }
    }
  });

  return true;
}

bool ChyokotovADenseMatMulFoxAlgorithmTBB::PostProcessingImpl() {
  return true;
}

}  // namespace chyokotov_a_dense_matrix_mul_foxs_algorithm
