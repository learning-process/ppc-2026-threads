#include "posternak_a_crs_mul_complex_matrix/seq/include/ops_seq.hpp"

#include <algorithm>
#include <complex>
#include <unordered_map>
#include <vector>

#include "posternak_a_crs_mul_complex_matrix/common/include/common.hpp"

namespace posternak_a_crs_mul_complex_matrix {

PosternakACRSMulComplexMatrixSEQ::PosternakACRSMulComplexMatrixSEQ(const InType &in) {
  SetTypeOfTask(GetStaticTypeOfTask());
  GetInput() = in;
  GetOutput() = CRSMatrix{};
}

bool PosternakACRSMulComplexMatrixSEQ::ValidationImpl() {
  const auto &[a, b] = GetInput();

  if (!a.IsValid() || !b.IsValid()) {
    return false;
  }

  if (a.cols != b.rows) {
    return false;
  }

  return true;
}

bool PosternakACRSMulComplexMatrixSEQ::PreProcessingImpl() {
  const auto &[a, b] = GetInput();

  GetOutput().rows = a.rows;
  GetOutput().cols = b.cols;
  GetOutput().values.clear();
  GetOutput().index_col.clear();
  GetOutput().index_row.clear();
  GetOutput().index_row.reserve(a.rows + 1);

  return true;
}

bool PosternakACRSMulComplexMatrixSEQ::RunImpl() {
  const auto &[a, b] = GetInput();
  auto &res = GetOutput();

  // одна из матриц пустая - результат пустой
  if (a.values.empty() || b.values.empty()) {
    res.index_row.assign(res.rows + 1, 0);
    return true;
  }

  for (int row = 0; row < res.rows; ++row) {
    // создаем хранилище результатов, где:
    // ключ - номер столбца в res
    // значение - сумма произведений столбца
    std::unordered_map<int, std::complex<double>> row_sum;

    for (int idx_a = a.index_row[row]; idx_a < a.index_row[row + 1]; ++idx_a) {
      const int col_a = a.index_col[idx_a];
      const auto &val_a = a.values[idx_a];

      for (int idx_b = b.index_row[col_a]; idx_b < b.index_row[col_a + 1]; ++idx_b) {
        const int col_b = b.index_col[idx_b];
        const auto &val_b = b.values[idx_b];

        row_sum[col_b] += val_a * val_b;
      }
    }

    // указатель на начало текущей строки в res
    res.index_row.push_back(static_cast<int>(res.values.size()));

    // копирование и сортировка результата
    std::vector<std::pair<int, std::complex<double>>> sorted_elements(row_sum.begin(), row_sum.end());
    std::sort(sorted_elements.begin(), sorted_elements.end(),
              [](const auto &a, const auto &b) { return a.first < b.first; });

    for (const auto &[col_idx, value] : sorted_elements) {
      // добавляем только значимые ненулевые элементы
      if (std::abs(value) > 1e-12) {
        res.values.push_back(value);
        res.index_col.push_back(col_idx);
      }
    }
  }

  // последнее значение index_row - общее количество элементов
  res.index_row.push_back(static_cast<int>(res.values.size()));

  return res.IsValid();
}

bool PosternakACRSMulComplexMatrixSEQ::PostProcessingImpl() {
  return GetOutput().IsValid();
}

}  // namespace posternak_a_crs_mul_complex_matrix
