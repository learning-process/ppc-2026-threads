#include "makoveeva_matmul_double_seq/seq/include/ops_seq.hpp"

#include <algorithm>
#include <cmath>
#include <vector>

namespace makoveeva_matmul_double_seq {
namespace {

void MultiplyBlocks(const std::vector<double> &a, const std::vector<double> &b, std::vector<double> &c, int n,
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

}  // namespace

MatmulDoubleSeqTask::MatmulDoubleSeqTask(const InType &in)
    : n_(std::get<0>(in)), A_(std::get<1>(in)), B_(std::get<2>(in)), C_(n_ * n_, 0.0) {
  SetTypeOfTask(GetStaticTypeOfTask());
  GetOutput() = C_;
}

bool MatmulDoubleSeqTask::ValidationImpl() {
  return n_ > 0 && A_.size() == n_ * n_ && B_.size() == n_ * n_;
}

bool MatmulDoubleSeqTask::PreProcessingImpl() {
  return true;
}

bool MatmulDoubleSeqTask::RunImpl() {
  if (n_ <= 0) {
    return false;
  }

  int n_int = static_cast<int>(n_);
  int block_size = std::max(1, static_cast<int>(std::sqrt(static_cast<double>(n_))));
  int num_blocks = (n_int + block_size - 1) / block_size;

  for (int ib = 0; ib < num_blocks; ++ib) {
    for (int jb = 0; jb < num_blocks; ++jb) {
      for (int kb = 0; kb < num_blocks; ++kb) {
        int i_start = ib * block_size;
        int i_end = std::min(i_start + block_size, n_int);
        int j_start = jb * block_size;
        int j_end = std::min(j_start + block_size, n_int);
        int k_start = kb * block_size;
        int k_end = std::min(k_start + block_size, n_int);

        MultiplyBlocks(A_, B_, C_, n_int, i_start, i_end, j_start, j_end, k_start, k_end);
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
