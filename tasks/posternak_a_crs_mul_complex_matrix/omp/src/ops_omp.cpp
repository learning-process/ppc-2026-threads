#include "posternak_a_crs_mul_complex_matrix/omp/include/ops_omp.hpp"

#include <algorithm>
#include <cmath>
#include <complex>
#include <cstddef>
#include <unordered_map>
#include <utility>
#include <vector>

#include "posternak_a_crs_mul_complex_matrix/common/include/common.hpp"

namespace {

size_t ComputeRowNoZeroCount(const posternak_a_crs_mul_complex_matrix::CRSMatrix &a,
                             const posternak_a_crs_mul_complex_matrix::CRSMatrix &b, int row, double threshold) {
  std::unordered_map<int, std::complex<double>> row_sum;

  for (int idx_a = a.index_row[row]; idx_a < a.index_row[row + 1]; ++idx_a) {
    int col_a = a.index_col[idx_a];
    auto val_a = a.values[idx_a];

    for (int idx_b = b.index_row[col_a]; idx_b < b.index_row[col_a + 1]; ++idx_b) {
      int col_b = b.index_col[idx_b];
      auto val_b = b.values[idx_b];
      row_sum[col_b] += val_a * val_b;
    }
  }

  size_t local = 0;
  for (const auto &[col, val] : row_sum) {
    if (std::abs(val) > threshold) {
      ++local;
    }
  }
  return local;
}

void BuildResultStructure(posternak_a_crs_mul_complex_matrix::CRSMatrix &res, std::vector<size_t> &row_prefix) {
  for (int i = 1; i < res.rows; ++i) {
    row_prefix[i] += row_prefix[i - 1];
  }

  const size_t total = row_prefix.empty() ? 0 : row_prefix.back();
  res.values.resize(total);
  res.index_col.resize(total);
  res.index_row.resize(res.rows + 1);

  for (int i = 0; i <= res.rows; ++i) {
    res.index_row[i] = (i == 0 ? 0 : static_cast<int>(row_prefix[i - 1]));
  }
}

void ComputeAndWriteRow(const posternak_a_crs_mul_complex_matrix::CRSMatrix &a,
                        const posternak_a_crs_mul_complex_matrix::CRSMatrix &b,
                        posternak_a_crs_mul_complex_matrix::CRSMatrix &res, int row, double threshold) {
  std::unordered_map<int, std::complex<double>> row_sum;

  for (int idx_a = a.index_row[row]; idx_a < a.index_row[row + 1]; ++idx_a) {
    int col_a = a.index_col[idx_a];
    auto val_a = a.values[idx_a];

    for (int idx_b = b.index_row[col_a]; idx_b < b.index_row[col_a + 1]; ++idx_b) {
      int col_b = b.index_col[idx_b];
      auto val_b = b.values[idx_b];
      row_sum[col_b] += val_a * val_b;
    }
  }

  std::vector<std::pair<int, std::complex<double>>> sorted(row_sum.begin(), row_sum.end());

  std::ranges::sort(sorted, [](const auto &p1, const auto &p2) { return p1.first < p2.first; });

  size_t pos = res.index_row[row];
  for (const auto &[col_idx, value] : sorted) {
    if (std::abs(value) > threshold) {
      res.values[pos] = value;
      res.index_col[pos] = col_idx;
      ++pos;
    }
  }
}

}  // namespace

namespace posternak_a_crs_mul_complex_matrix {

PosternakACRSMulComplexMatrixOMP::PosternakACRSMulComplexMatrixOMP(const InType &in) {
  SetTypeOfTask(GetStaticTypeOfTask());
  GetInput() = in;
  GetOutput() = CRSMatrix{};
}

bool PosternakACRSMulComplexMatrixOMP::ValidationImpl() {
  const auto &input = GetInput();
  const auto &a = input.first;
  const auto &b = input.second;
  return a.IsValid() && b.IsValid() && a.cols == b.rows;
}

bool PosternakACRSMulComplexMatrixOMP::PreProcessingImpl() {
  const auto &input = GetInput();
  const auto &a = input.first;
  const auto &b = input.second;
  auto &res = GetOutput();

  res.rows = a.rows;
  res.cols = b.cols;
  return true;
}

bool PosternakACRSMulComplexMatrixOMP::RunImpl() {
  const auto &input = GetInput();
  const auto &a = input.first;
  const auto &b = input.second;
  auto &res = GetOutput();

  if (a.values.empty() || b.values.empty()) {
    res.values.clear();
    res.index_col.clear();
    res.index_row.assign(res.rows + 1, 0);
    return true;
  }

  constexpr double kThreshold = 1e-12;
  std::vector<size_t> no_zero_rows(res.rows, 0);

  // каждый поток определяет количество ненулевых элементов в своих строках
#pragma omp parallel for default(none) shared(a, b, res, no_zero_rows) schedule(dynamic)
  for (int row = 0; row < res.rows; ++row) {
    no_zero_rows[row] = ComputeRowNoZeroCount(a, b, row, kThreshold);
  }

  // структурируем результат, чтобы избежать конфликта потоков (предотвращаем гонку данных)
#pragma omp single
  {
    BuildResultStructure(res, no_zero_rows);
  }

  // записываем результат в итоговый массив параллельно
#pragma omp parallel for default(none) shared(a, b, res) schedule(dynamic)
  for (int row = 0; row < res.rows; ++row) {
    ComputeAndWriteRow(a, b, res, row, kThreshold);
  }

  return res.IsValid();
}

bool PosternakACRSMulComplexMatrixOMP::PostProcessingImpl() {
  return GetOutput().IsValid();
}

}  // namespace posternak_a_crs_mul_complex_matrix
