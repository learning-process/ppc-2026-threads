#include "makoveeva_matmul_double_seq/seq/include/ops_seq.hpp"

#include <algorithm>
#include <cmath>
#include <vector>

#include "makoveeva_matmul_double_seq/common/include/common.hpp"  // для InType

namespace makoveeva_matmul_double_seq {
namespace {

void ProcessBlock(const std::vector<double> &a, const std::vector<double> &b, std::vector<double> &c, int n,
                  int i_start, int i_end, int j_start, int j_end, int k_start, int k_end) {
  for (int i = i_start; i < i_end; ++i) {
    for (int j = j_start; j < j_end; ++j) {
      double sum = 0.0;
      for (int k = k_start; k < k_end; ++k) {
        sum += a[(i * n) + k] * b[(k * n) + j];
      }
      c[(i * n) + j] += sum;
    }
  }
}

int CalculateBlockSize(int n) {
  return std::max(1, static_cast<int>(std::sqrt(static_cast<double>(n))));
}

int CalculateNumBlocks(int n, int block_size) {
  return (n + block_size - 1) / block_size;
}

}  // namespace

MatmulDoubleSeqTask::MatmulDoubleSeqTask(const InType &in)
    : n_(std::get<0>(in)), A_(std::get<1>(in)), B_(std::get<2>(in)), C_(n_ * n_, 0.0) {
  SetTypeOfTask(GetStaticTypeOfTask());
  GetOutput() = C_;
}

bool MatmulDoubleSeqTask::ValidationImpl() {
  const bool valid_n = n_ > 0;
  const bool valid_a = A_.size() == n_ * n_;
  const bool valid_b = B_.size() == n_ * n_;
  return valid_n && valid_a && valid_b;
}

bool MatmulDoubleSeqTask::PreProcessingImpl() {
  return true;
}

bool MatmulDoubleSeqTask::RunImpl() {
  if (n_ <= 0) {
    return false;
  }

  // Очищаем C_ перед вычислениями
  C_.assign(C_.size(), 0.0);

  const int n_int = static_cast<int>(n_);
  const int block_size = CalculateBlockSize(n_int);
  const int num_blocks = CalculateNumBlocks(n_int, block_size);

  for (int ib = 0; ib < num_blocks; ++ib) {
    for (int jb = 0; jb < num_blocks; ++jb) {
      for (int kb = 0; kb < num_blocks; ++kb) {
        const int i_start = ib * block_size;
        const int i_end = std::min(i_start + block_size, n_int);
        const int j_start = jb * block_size;
        const int j_end = std::min(j_start + block_size, n_int);
        const int k_start = kb * block_size;
        const int k_end = std::min(k_start + block_size, n_int);

        ProcessBlock(A_, B_, C_, n_int, i_start, i_end, j_start, j_end, k_start, k_end);
      }
    }
  }

  GetOutput() = C_;
  return true;
}

bool MatmulDoubleSeqTask::PostProcessingImpl() {
  return true;
}

}  // namespace makoveeva_matmul_double_seq
