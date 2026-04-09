#include "volkov_a_sparse_mat_mul_ccs/tbb/include/ops_tbb.hpp"

#include <tbb/blocked_range.h>
#include <tbb/parallel_for.h>

#include <cmath>
#include <tuple>
#include <vector>

#include "volkov_a_sparse_mat_mul_ccs/common/include/common.hpp"

namespace volkov_a_sparse_mat_mul_ccs {

VolkovASparseMatMulCcsTbb::VolkovASparseMatMulCcsTbb(const InType &in) {
  SetTypeOfTask(GetStaticTypeOfTask());
  GetInput() = in;
}

bool VolkovASparseMatMulCcsTbb::ValidationImpl() {
  // Используем вашу логику валидации из SEQ версии
  const auto &matrix_a = std::get<0>(GetInput());
  const auto &matrix_b = std::get<1>(GetInput());
  return (matrix_a.cols_count == matrix_b.rows_count);
}

bool VolkovASparseMatMulCcsTbb::PreProcessingImpl() {
  return true;
}

bool VolkovASparseMatMulCcsTbb::RunImpl() {
  const auto &matrix_a = std::get<0>(GetInput());
  const auto &matrix_b = std::get<1>(GetInput());
  auto &matrix_c = GetOutput();

  matrix_c.rows_count = matrix_a.rows_count;
  matrix_c.cols_count = matrix_b.cols_count;

  // Локальные хранилища для каждого столбца, чтобы избежать гонок данных (data races)
  std::vector<std::vector<int>> local_row_indices(matrix_b.cols_count);
  std::vector<std::vector<double>> local_values(matrix_b.cols_count);

  // ПЕРВЫЙ ЭТАП: Параллельное вычисление значений столбцов
  tbb::parallel_for(
      tbb::blocked_range<int>(0, matrix_b.cols_count), [&](const tbb::blocked_range<int> &range) {
        // Локальный аккумулятор вынесен внутрь потока
        std::vector<double> col_accumulator(matrix_a.rows_count, 0.0);

        for (int j = range.begin(); j != range.end(); ++j) {
          int b_start = matrix_b.col_ptrs[j];
          int b_end = matrix_b.col_ptrs[j + 1];

          for (int k = b_start; k < b_end; ++k) {
            int b_row = matrix_b.row_indices[k];
            double b_val = matrix_b.values[k];

            int a_start = matrix_a.col_ptrs[b_row];
            int a_end = matrix_a.col_ptrs[b_row + 1];

            for (int idx = a_start; idx < a_end; ++idx) {
              int a_row = matrix_a.row_indices[idx];
              double a_val = matrix_a.values[idx];
              col_accumulator[a_row] += a_val * b_val;
            }
          }

          // Сбор ненулевых элементов для текущего столбца j
          for (int i = 0; i < matrix_a.rows_count; ++i) {
            if (std::abs(col_accumulator[i]) > 1e-10) {
              local_row_indices[j].push_back(i);
              local_values[j].push_back(col_accumulator[i]);
            }
            col_accumulator[i] = 0.0;  // Сброс аккумулятора для следующей итерации
          }
        }
      });

  // ВТОРОЙ ЭТАП: Последовательный подсчет указателей столбцов (префиксная сумма)
  matrix_c.col_ptrs.assign(matrix_c.cols_count + 1, 0);
  for (int j = 0; j < matrix_b.cols_count; ++j) {
    matrix_c.col_ptrs[j + 1] = matrix_c.col_ptrs[j] + static_cast<int>(local_row_indices[j].size());
  }

  // Выделение памяти под итоговые массивы
  matrix_c.non_zeros = matrix_c.col_ptrs.back();
  matrix_c.row_indices.resize(matrix_c.non_zeros);
  matrix_c.values.resize(matrix_c.non_zeros);

  // ТРЕТИЙ ЭТАП: Параллельное копирование данных из локальных векторов в глобальные
  tbb::parallel_for(tbb::blocked_range<int>(0, matrix_b.cols_count),
                    [&](const tbb::blocked_range<int> &range) {
                      for (int j = range.begin(); j != range.end(); ++j) {
                        int offset = matrix_c.col_ptrs[j];
                        int current_col_size = static_cast<int>(local_row_indices[j].size());

                        for (int k = 0; k < current_col_size; ++k) {
                          matrix_c.row_indices[offset + k] = local_row_indices[j][k];
                          matrix_c.values[offset + k] = local_values[j][k];
                        }
                      }
                    });

  return true;
}

bool VolkovASparseMatMulCcsTbb::PostProcessingImpl() {
  return true;
}

}  // namespace volkov_a_sparse_mat_mul_ccs