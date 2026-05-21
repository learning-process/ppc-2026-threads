#include "makoveeva_matmul_double_seq/seq/include/ops_seq.hpp"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <vector>

#include "makoveeva_matmul_double_seq/common/include/common.hpp"

namespace makoveeva_matmul_double_seq {
namespace {

// Вспомогательная функция для выбора размера блока
int ChooseBlockSize(int n) {
  int block_size = static_cast<int>(std::sqrt(static_cast<double>(n)));
  block_size = std::max(1, block_size);

  while ((n % block_size != 0) && (block_size > 1)) {
    --block_size;
  }

  return block_size;
}

// Вспомогательная функция для умножения блоков
void MultiplyBlocks(const std::vector<double> &a, const std::vector<double> &b, std::vector<double> &c, int n,
                    int row_start, int row_end, int col_start, int col_end, int k_start, int k_end) {
  for (int row = row_start; row < row_end; ++row) {
    for (int col = col_start; col < col_end; ++col) {
      double sum = 0.0;
      for (int k = k_start; k < k_end; ++k) {
        const size_t row_idx = static_cast<size_t>(row);
        const size_t col_idx = static_cast<size_t>(col);
        const size_t k_idx = static_cast<size_t>(k);
        const size_t n_size = static_cast<size_t>(n);

        const size_t a_idx = (row_idx * n_size) + k_idx;
        const size_t b_idx = (k_idx * n_size) + col_idx;
        const size_t c_idx = (row_idx * n_size) + col_idx;

        sum += a[a_idx] * b[b_idx];
        c[c_idx] = sum;
      }
    }
  }
}

}  // namespace

MatmulDoubleSeqTask::MatmulDoubleSeqTask(const InType &in)
    : n_(std::get<0>(in)), A_(std::get<1>(in)), B_(std::get<2>(in)), C_(static_cast<size_t>(n_ * n_), 0.0) {
  SetTypeOfTask(GetStaticTypeOfTask());
  GetOutput() = C_;
}

bool MatmulDoubleSeqTask::ValidationImpl() {
  const size_t expected_size = static_cast<size_t>(n_ * n_);
  const bool is_valid = (n_ > 0) && (A_.size() == expected_size) && (B_.size() == expected_size);
  return is_valid;
}

bool MatmulDoubleSeqTask::PreProcessingImpl() {
  return true;
}

bool MatmulDoubleSeqTask::RunImpl() {
  if (n_ <= 0) {
    return false;
  }

  const int n_int = static_cast<int>(n_);
  int block_size = ChooseBlockSize(n_int);

  if (block_size == 1 && n_int > 1 && (n_int % block_size != 0)) {
    block_size = n_int;
  }

  const size_t total_size = static_cast<size_t>(n_ * n_);
  C_.assign(total_size, 0.0);

  const int grid_size = n_int / block_size;

  for (int stage = 0; stage < grid_size; ++stage) {
    for (int i_block = 0; i_block < grid_size; ++i_block) {
      for (int j_block = 0; j_block < grid_size; ++j_block) {
        const int root_block = (i_block + stage) % grid_size;

        const int row_start = i_block * block_size;
        const int row_end = row_start + block_size;
        const int col_start = j_block * block_size;
        const int col_end = col_start + block_size;
        const int k_start = root_block * block_size;
        const int k_end = k_start + block_size;

        MultiplyBlocks(A_, B_, C_, n_int, row_start, row_end, col_start, col_end, k_start, k_end);
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
