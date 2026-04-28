#include "volkov_a_sparse_mat_mul_ccs/omp/include/ops_omp.hpp"

#include <omp.h>

#include <cmath>
#include <tuple>
#include <vector>

#include "volkov_a_sparse_mat_mul_ccs/common/include/common.hpp"

namespace volkov_a_sparse_mat_mul_ccs {

VolkovASparseMatMulCcsOmp::VolkovASparseMatMulCcsOmp(const InType &in) {
  SetTypeOfTask(GetStaticTypeOfTask());
  GetInput() = in;
}

bool VolkovASparseMatMulCcsOmp::ValidationImpl() {
  const auto &matrix_a = std::get<0>(GetInput());
  const auto &matrix_b = std::get<1>(GetInput());
  return (matrix_a.cols_count == matrix_b.rows_count);
}

bool VolkovASparseMatMulCcsOmp::PreProcessingImpl() {
  return true;
}

bool VolkovASparseMatMulCcsOmp::RunImpl() {
  const auto &matrix_a = std::get<0>(GetInput());
  const auto &matrix_b = std::get<1>(GetInput());
  auto &matrix_c = GetOutput();

  matrix_c.rows_count = matrix_a.rows_count;
  matrix_c.cols_count = matrix_b.cols_count;
  matrix_c.col_ptrs.assign(matrix_c.cols_count + 1, 0);

  // Временные массивы для хранения результатов каждого столбца до их слияния
  std::vector<std::vector<int>> cols_res_rows(matrix_b.cols_count);
  std::vector<std::vector<double>> cols_res_values(matrix_b.cols_count);

  // Параллельный регион. Всё, что объявлено внутри, становится локальным для потока!
#pragma omp parallel default(none) shared(matrix_a, matrix_b, cols_res_rows, cols_res_values)
  {
    // Изящный подход: потоколокальный аккумулятор.
    // Никаких omp_get_thread_num() и сложной индексации!
    std::vector<double> local_accumulator(matrix_a.rows_count, 0.0);

    // Распределяем итерации по потокам (dynamic лучше для разреженных матриц)
#pragma omp for schedule(dynamic)
    for (int j = 0; j < matrix_b.cols_count; ++j) {
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

          local_accumulator[a_row] += a_val * b_val;
        }
      }

      // Переписываем ненулевые элементы из локального аккумулятора в вектор столбца
      for (int i = 0; i < matrix_a.rows_count; ++i) {
        if (std::abs(local_accumulator[i]) > 1e-10) {
          cols_res_rows[j].push_back(i);
          cols_res_values[j].push_back(local_accumulator[i]);
        }
        // Сразу обнуляем для вычисления следующего столбца в этом же потоке
        local_accumulator[i] = 0.0;
      }
    }
  }

  // Однопоточное вычисление префиксных сумм для col_ptrs
  int current_non_zeros = 0;
  for (int j = 0; j < matrix_b.cols_count; ++j) {
    matrix_c.col_ptrs[j] = current_non_zeros;
    current_non_zeros += static_cast<int>(cols_res_values[j].size());
  }
  matrix_c.col_ptrs[matrix_b.cols_count] = current_non_zeros;
  matrix_c.non_zeros = current_non_zeros;

  matrix_c.row_indices.resize(current_non_zeros);
  matrix_c.values.resize(current_non_zeros);

  // Параллельное копирование временных данных в итоговый одномерный массив формата CCS
#pragma omp parallel for default(none) schedule(static) shared(matrix_b, matrix_c, cols_res_rows, cols_res_values)
  for (int j = 0; j < matrix_b.cols_count; ++j) {
    int offset = matrix_c.col_ptrs[j];
    int elements_in_col = static_cast<int>(cols_res_values[j].size());
    for (int k = 0; k < elements_in_col; ++k) {
      matrix_c.row_indices[offset + k] = cols_res_rows[j][k];
      matrix_c.values[offset + k] = cols_res_values[j][k];
    }
  }

  return true;
}

bool VolkovASparseMatMulCcsOmp::PostProcessingImpl() {
  return true;
}

}  // namespace volkov_a_sparse_mat_mul_ccs
