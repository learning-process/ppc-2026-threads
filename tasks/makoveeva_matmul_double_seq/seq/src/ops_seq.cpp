#include "makoveeva_matmul_double_seq/seq/include/ops_seq.hpp"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <vector>

#include "makoveeva_matmul_double_seq/common/include/common.hpp"

namespace makoveeva_matmul_double_seq {

MatmulDoubleSeqTask::MatmulDoubleSeqTask(const InType &in)
    : n_(std::get<0>(in)), A_(std::get<1>(in)), B_(std::get<2>(in)), C_(static_cast<std::size_t>(n_ * n_), 0.0) {
  SetTypeOfTask(GetStaticTypeOfTask());
  GetOutput() = C_;
}

bool MatmulDoubleSeqTask::ValidationImpl() {
  const std::size_t expected_size = static_cast<std::size_t>(n_ * n_);
  const bool is_valid = n_ > 0 && A_.size() == expected_size && B_.size() == expected_size;
  return is_valid;
}

bool MatmulDoubleSeqTask::PreProcessingImpl() {
  return true;
}

namespace {

// Helper function to get linear index
std::size_t Idx(int i, int j, int dim) {
  return static_cast<std::size_t>(i * dim + j);
}

// Helper function to multiply block A_block * B_block and add to C_block
void MultiplyAndAddBlock(const std::vector<double> &A_block, const std::vector<double> &B_block,
                         std::vector<double> &C_block, int block_size) {
  for (int i = 0; i < block_size; ++i) {
    for (int k = 0; k < block_size; ++k) {
      const double a_ik = A_block[Idx(i, k, block_size)];
      for (int j = 0; j < block_size; ++j) {
        C_block[Idx(i, j, block_size)] += a_ik * B_block[Idx(k, j, block_size)];
      }
    }
  }
}

// Helper function to copy a block from matrix
std::vector<double> GetBlock(const std::vector<double> &matrix, int block_row, int block_col, int block_size,
                             int matrix_dim) {
  std::vector<double> block(static_cast<std::size_t>(block_size * block_size), 0.0);
  for (int i = 0; i < block_size; ++i) {
    const int global_i = block_row * block_size + i;
    if (global_i >= matrix_dim) {
      continue;
    }

    for (int j = 0; j < block_size; ++j) {
      const int global_j = block_col * block_size + j;
      if (global_j >= matrix_dim) {
        continue;
      }

      block[Idx(i, j, block_size)] = matrix[Idx(global_i, global_j, matrix_dim)];
    }
  }
  return block;
}

// Helper function to add a block to matrix
void AddBlockToMatrix(const std::vector<double> &block, std::vector<double> &matrix, int block_row, int block_col,
                      int block_size, int matrix_dim) {
  for (int i = 0; i < block_size; ++i) {
    const int global_i = block_row * block_size + i;
    if (global_i >= matrix_dim) {
      continue;
    }

    for (int j = 0; j < block_size; ++j) {
      const int global_j = block_col * block_size + j;
      if (global_j >= matrix_dim) {
        continue;
      }

      matrix[Idx(global_i, global_j, matrix_dim)] += block[Idx(i, j, block_size)];
    }
  }
}

// Fox's algorithm core function
void FoxAlgorithm(const std::vector<double> &A, const std::vector<double> &B, std::vector<double> &C, int matrix_dim,
                  int block_size) {
  const int grid_size = matrix_dim / block_size;

  // Initialize result matrix C with zeros
  std::fill(C.begin(), C.end(), 0.0);

  // Main loop of Fox's algorithm
  for (int stage = 0; stage < grid_size; ++stage) {
    // For each row block
    for (int row_block = 0; row_block < grid_size; ++row_block) {
      // The diagonal block of A that will be broadcast in this stage
      const int a_block_col = (row_block + stage) % grid_size;

      // Get the block from A
      std::vector<double> A_block = GetBlock(A, row_block, a_block_col, block_size, matrix_dim);

      // For each column block
      for (int col_block = 0; col_block < grid_size; ++col_block) {
        // Get the block from B
        std::vector<double> B_block = GetBlock(B, a_block_col, col_block, block_size, matrix_dim);

        // Get current C block
        std::vector<double> C_block = GetBlock(C, row_block, col_block, block_size, matrix_dim);

        // Multiply A_block * B_block and add to C_block
        MultiplyAndAddBlock(A_block, B_block, C_block, block_size);

        // Put C_block back
        AddBlockToMatrix(C_block, C, row_block, col_block, block_size, matrix_dim);
      }
    }
  }
}

}  // namespace

bool MatmulDoubleSeqTask::RunImpl() {
  if (n_ <= 0) {
    return false;
  }

  const int matrix_dim = static_cast<int>(n_);

  // Calculate optimal block size for Fox's algorithm
  int block_size = static_cast<int>(std::sqrt(static_cast<double>(matrix_dim)));
  block_size = std::max(1, block_size);

  // Ensure block_size divides matrix_dim evenly
  while (matrix_dim % block_size != 0 && block_size > 1) {
    --block_size;
  }

  // If no suitable block size found, use matrix_dim itself as block size
  if (block_size == 1 && matrix_dim > 1 && matrix_dim % block_size != 0) {
    block_size = matrix_dim;
  }

  // Allocate result matrix
  C_.assign(static_cast<std::size_t>(matrix_dim * matrix_dim), 0.0);

  // Execute Fox's algorithm
  FoxAlgorithm(A_, B_, C_, matrix_dim, block_size);

  GetOutput() = C_;
  return true;
}

bool MatmulDoubleSeqTask::PostProcessingImpl() {
  return true;
}

}  // namespace makoveeva_matmul_double_seq
