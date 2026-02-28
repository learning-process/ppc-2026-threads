#include "chyokotov_a_dense_matrix_mul_foxs_algorithm/seq/include/ops_seq.hpp"

#include <algorithm>
#include <cmath>
#include <vector>

#include "chyokotov_a_dense_matrix_mul_foxs_algorithm/common/include/common.hpp"

namespace chyokotov_a_dense_matrix_mul_foxs_algorithm {

ChyokotovADenseMatMulFoxAlgorithmSEQ::ChyokotovADenseMatMulFoxAlgorithmSEQ(const InType &in) {
  SetTypeOfTask(GetStaticTypeOfTask());
  GetInput() = in;
  GetOutput().clear();
}

bool ChyokotovADenseMatMulFoxAlgorithmSEQ::ValidationImpl() {
  return (GetInput().first.size() == GetInput().second.size());
}

bool ChyokotovADenseMatMulFoxAlgorithmSEQ::PreProcessingImpl() {
  GetOutput().clear();
  GetOutput().resize(GetInput().first.size(), 0.0);
  return true;
}

int ChyokotovADenseMatMulFoxAlgorithmSEQ::CalculateBlockSize(int n) {
  return static_cast<int>(std::sqrt(static_cast<double>(n)));
}

int ChyokotovADenseMatMulFoxAlgorithmSEQ::CountBlock(int n, int size) {
  return (n + size - 1) / size;
}

void ChyokotovADenseMatMulFoxAlgorithmSEQ::Matmul(std::vector<double> &a, std::vector<double> &b, int n, int istart,
                                                  int iend, int jstart, int jend, int kstart, int kend) {
  for (int i = istart; i < iend; i++) {
    for (int j = jstart; j < jend; j++) {
      double sum = 0.0;
      for (int k = kstart; k < kend; k++) {
        sum += a[(i * n) + k] * b[(k * n) + j];
      }
      GetOutput()[(i * n) + j] += sum;
    }
  }
}

bool ChyokotovADenseMatMulFoxAlgorithmSEQ::RunImpl() {
  std::vector<double> a = GetInput().first;
  std::vector<double> b = GetInput().second;
  int n = static_cast<int>(std::sqrt(static_cast<double>(a.size())));
  if (n == 0) {
    return true;
  }

  int block_size = CalculateBlockSize(n);
  int count_block = CountBlock(n, block_size);

  for (int ic = 0; ic < count_block; ic++) {
    for (int jc = 0; jc < count_block; jc++) {
      for (int kc = 0; kc < count_block; kc++) {
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

  return true;
}

bool ChyokotovADenseMatMulFoxAlgorithmSEQ::PostProcessingImpl() {
  return true;
}

}  // namespace chyokotov_a_dense_matrix_mul_foxs_algorithm
