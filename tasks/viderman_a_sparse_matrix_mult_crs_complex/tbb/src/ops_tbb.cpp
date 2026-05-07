#include "viderman_a_sparse_matrix_mult_crs_complex/tbb/include/ops_tbb.hpp"

#include <algorithm>
#include <cmath>
#include <vector>

#include "oneapi/tbb/enumerable_thread_specific.h"
#include "oneapi/tbb/parallel_for.h"

namespace viderman_a_sparse_matrix_mult_crs_complex {

// Структура для локальных буферов каждого потока
struct ThreadLocalBuffers {
  std::vector<Complex> accumulator;
  std::vector<int> marker;
  std::vector<int> active_indices;

  explicit ThreadLocalBuffers(int size) : accumulator(size, Complex(0.0, 0.0)), marker(size, -1) {}
};

VidermanASparseMatrixMultCRSComplexTBB::VidermanASparseMatrixMultCRSComplexTBB(const InType &in) {
  SetTypeOfTask(GetStaticTypeOfTask());
  GetInput() = in;
  // Инициализируем выходную матрицу (OutType здесь это CRSMatrix)
  GetOutput() = CRSMatrix();
}

bool VidermanASparseMatrixMultCRSComplexTBB::ValidationImpl() {
  const auto &input = GetInput();
  const auto &a = std::get<0>(input);
  const auto &b = std::get<1>(input);
  return a.IsValid() && b.IsValid() && (a.cols == b.rows);
}

bool VidermanASparseMatrixMultCRSComplexTBB::PreProcessingImpl() {
  const auto &input = GetInput();
  a_ = &std::get<0>(input);
  b_ = &std::get<1>(input);
  return true;
}

bool VidermanASparseMatrixMultCRSComplexTBB::RunImpl() {
  // ИСПРАВЛЕНО: Передаем GetOutput() напрямую без std::get
  Multiply(*a_, *b_, GetOutput());
  return true;
}

bool VidermanASparseMatrixMultCRSComplexTBB::PostProcessingImpl() {
  return true;
}

void VidermanASparseMatrixMultCRSComplexTBB::Multiply(const CRSMatrix &a, const CRSMatrix &b, CRSMatrix &c) {
  c.rows = a.rows;
  c.cols = b.cols;
  c.row_ptr.assign(a.rows + 1, 0);

  // Временные контейнеры для строк
  std::vector<std::vector<int>> all_col_indices(a.rows);
  std::vector<std::vector<Complex>> all_values(a.rows);

  // TLS для исключения Race Condition и лишних аллокаций
  oneapi::tbb::enumerable_thread_specific<ThreadLocalBuffers> tls([&]() { return ThreadLocalBuffers(b.cols); });

  oneapi::tbb::parallel_for(0, a.rows, [&](int i) {
    auto &loc = tls.local();
    auto &acc = loc.accumulator;
    auto &mark = loc.marker;
    auto &active = loc.active_indices;

    active.clear();

    for (int j = a.row_ptr[i]; j < a.row_ptr[i + 1]; ++j) {
      int col_a = a.col_indices[j];
      Complex val_a = a.values[j];

      for (int k = b.row_ptr[col_a]; k < b.row_ptr[col_a + 1]; ++k) {
        int col_b = b.col_indices[k];
        acc[col_b] += val_a * b.values[k];

        if (mark[col_b] != i) {
          mark[col_b] = i;
          active.push_back(col_b);
        }
      }
    }

    // Сортировка для CRS формата
    std::sort(active.begin(), active.end());

    for (int col : active) {
      if (std::abs(acc[col]) > 1e-15) {  // Используем порог из ваших тестов
        all_col_indices[i].push_back(col);
        all_values[i].push_back(acc[col]);
      }
      acc[col] = Complex(0.0, 0.0);  // Очистка
    }
  });

  // Сборка результирующей матрицы
  for (int i = 0; i < a.rows; ++i) {
    c.row_ptr[i + 1] = c.row_ptr[i] + static_cast<int>(all_col_indices[i].size());
  }

  c.col_indices.reserve(c.row_ptr[a.rows]);
  c.values.reserve(c.row_ptr[a.rows]);

  for (int i = 0; i < a.rows; ++i) {
    c.col_indices.insert(c.col_indices.end(), all_col_indices[i].begin(), all_col_indices[i].end());
    c.values.insert(c.values.end(), all_values[i].begin(), all_values[i].end());
  }
}

}  // namespace viderman_a_sparse_matrix_mult_crs_complex
