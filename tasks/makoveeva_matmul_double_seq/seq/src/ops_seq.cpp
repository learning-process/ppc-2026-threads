#include "makoveeva_matmul_double_seq/seq/include/ops_seq.hpp"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <vector>

#include "makoveeva_matmul_double_seq/common/include/common.hpp"

namespace makoveeva_matmul_double_seq {

MatmulDoubleSeqTask::MatmulDoubleSeqTask(const InType &in)
    : n_(std::get<0>(in)), A_(std::get<1>(in)), B_(std::get<2>(in)), C_(n_ * n_, 0.0) {
  SetTypeOfTask(GetStaticTypeOfTask());
  GetOutput() = C_;
}

bool MatmulDoubleSeqTask::ValidationImpl() {
  // Remove redundant casting - n_ * n_ already produces size_t
  const bool is_valid =
      n_ > 0 && A_.size() == static_cast<std::size_t>(n_ * n_) && B_.size() == static_cast<std::size_t>(n_ * n_);
  return is_valid;
}

bool MatmulDoubleSeqTask::PreProcessingImpl() {
  return true;
}

// Helper function to reduce cognitive complexity
namespace {
void ComputeBlock(const std::vector<double> &A, const std::vector<double> &B, std::vector<double> &C, int matrix_dim,
                  int row_start, int col_start, int root_start, int block_size) {
  for (int local_i = 0; local_i < block_size; ++local_i) {
    const int global_i = row_start + local_i;
    if (global_i >= matrix_dim) {
      continue;
    }

    for (int local_j = 0; local_j < block_size; ++local_j) {
      const int global_j = col_start + local_j;
      if (global_j >= matrix_dim) {
        continue;
      }

      double accumulator = 0.0;
      for (int local_k = 0; local_k < block_size; ++local_k) {
        const int global_k = root_start + local_k;
        if (global_k >= matrix_dim) {
          continue;
        }

        accumulator += A[static_cast<std::size_t>(global_i * matrix_dim + global_k)] *
                       B[static_cast<std::size_t>(global_k * matrix_dim + global_j)];
      }
      C[static_cast<std::size_t>(global_i * matrix_dim + global_j)] += accumulator;
    }
  }
}
}  // namespace

bool MatmulDoubleSeqTask::RunImpl() {
  if (n_ <= 0) {
    return false;
  }

  C_.assign(C_.size(), 0.0);

  const int matrix_dim = static_cast<int>(n_);
  int block_size = static_cast<int>(std::sqrt(static_cast<double>(matrix_dim)));
  block_size = std::max(1, block_size);

  // Adjust block size to divide matrix_dim evenly
  while (matrix_dim % block_size != 0 && block_size > 1) {
    --block_size;
  }

  const int grid_size = matrix_dim / block_size;

  // Matrix multiplication with blocking optimization
  for (int stage = 0; stage < grid_size; ++stage) {
    for (int row_block = 0; row_block < grid_size; ++row_block) {
      const int root_block = (row_block + stage) % grid_size;

      for (int col_block = 0; col_block < grid_size; ++col_block) {
        const int row_start = row_block * block_size;
        const int col_start = col_block * block_size;
        const int root_start = root_block * block_size;

        ComputeBlock(A_, B_, C_, matrix_dim, row_start, col_start, root_start, block_size);
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
