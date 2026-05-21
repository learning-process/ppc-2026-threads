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
  const auto expected_size = static_cast<std::size_t>(n_ * n_);
  const bool is_valid = n_ > 0 && A_.size() == expected_size && B_.size() == expected_size;
  return is_valid;
}

bool MatmulDoubleSeqTask::PreProcessingImpl() {
  return true;
}

bool MatmulDoubleSeqTask::RunImpl() {
  if (n_ <= 0) {
    return false;
  }

  // Блокировка для алгоритма Фокса
  int block_size = static_cast<int>(std::sqrt(static_cast<double>(n_)));
  block_size = std::max(1, block_size);
  
  // Убеждаемся, что размер матрицы делится на block_size
  while (n_ % block_size != 0 && block_size > 1) {
    --block_size;
  }
  
  if (block_size == 1 && n_ > 1 && n_ % block_size != 0) {
    block_size = n_;
  }

  // Инициализируем результат нулями
  C_.assign(static_cast<std::size_t>(n_ * n_), 0.0);
  
  const int grid_size = n_ / block_size;
  
  // Алгоритм Фокса
  for (int stage = 0; stage < grid_size; ++stage) {
    for (int i = 0; i < grid_size; ++i) {
      for (int j = 0; j < grid_size; ++j) {
        // Индекс диагонального блока A на этом этапе
        int root = (i + stage) % grid_size;
        
        // Перемножаем блоки
        for (int bi = 0; bi < block_size; ++bi) {
          for (int bj = 0; bj < block_size; ++bj) {
            for (int bk = 0; bk < block_size; ++bk) {
              int row = i * block_size + bi;
              int col = j * block_size + bj;
              int k_idx = root * block_size + bk;
              
              C_[row * n_ + col] += A_[row * n_ + k_idx] * B_[k_idx * n_ + col];
            }
          }
        }
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