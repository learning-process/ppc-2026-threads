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
  const auto n_size = static_cast<size_t>(n);

  for (int row = row_start; row < row_end; ++row) {
    const auto row_idx = static_cast<size_t>(row);
    const auto row_offset = row_idx * n_size;

    for (int col = col_start; col < col_end; ++col) {
      const auto col_idx = static_cast<size_t>(col);
      double sum = 0.0;

      for (int k = k_start; k < k_end; ++k) {
        const auto k_idx = static_cast<size_t>(k);
        const auto a_idx = row_offset + k_idx;
        const auto b_idx = (k_idx * n_size) + col_idx;
        sum += a[a_idx] * b[b_idx];
      }

      const auto c_idx = row_offset + col_idx;
      c[c_idx] += sum;
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
  const auto expected_size = n_ * n_;
  const auto is_valid = (n_ > 0) && (A_.size() == expected_size) && (B_.size() == expected_size);
  return is_valid;
}

bool MatmulDoubleSeqTask::PreProcessingImpl() {
  return true;
}

bool MatmulDoubleSeqTask::RunImpl() {
  if (n_ <= 0) {
    return false;
  }

  const auto n_int = static_cast<int>(n_);
  auto block_size = ChooseBlockSize(n_int);

  if (block_size == 1 && n_int > 1 && (n_int % block_size != 0)) {
    block_size = n_int;
  }

  const auto total_size = n_ * n_;
  C_.assign(total_size, 0.0);

  const auto grid_size = n_int / block_size;

  for (int stage = 0; stage < grid_size; ++stage) {
    for (int i_block = 0; i_block < grid_size; ++i_block) {
      for (int j_block = 0; j_block < grid_size; ++j_block) {
        const auto root_block = (i_block + stage) % grid_size;

        const auto row_start = i_block * block_size;
        const auto row_end = row_start + block_size;
        const auto col_start = j_block * block_size;
        const auto col_end = col_start + block_size;
        const auto k_start = root_block * block_size;
        const auto k_end = k_start + block_size;

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
