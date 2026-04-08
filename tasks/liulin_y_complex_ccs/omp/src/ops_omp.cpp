#include "liulin_y_complex_ccs/omp/include/ops_omp.hpp"

#include <omp.h>

#include <algorithm>
#include <cmath>
#include <complex>
#include <cstddef>
#include <vector>

#include "liulin_y_complex_ccs/common/include/common.hpp"

namespace liulin_y_complex_ccs {

namespace {
constexpr double kZeroThreshold = 1e-15;

bool IsValueNonZero(const std::complex<double> &value) {
  return std::abs(value) > kZeroThreshold;
}

void ProcessSingleColumn(int col_idx_b, const CCSMatrix &mat_a, const CCSMatrix &mat_b,
                         std::vector<std::complex<double>> &accumulator, std::vector<int> &active_rows,
                         std::vector<int> &row_marker, std::vector<std::complex<double>> &out_values,
                         std::vector<int> &out_rows) {
  const int col_start = mat_b.col_index[static_cast<size_t>(col_idx_b)];
  const int col_end = mat_b.col_index[static_cast<size_t>(col_idx_b) + 1];

  for (int idx_b = col_start; idx_b < col_end; ++idx_b) {
    const int row_idx_b = mat_b.row_index[static_cast<size_t>(idx_b)];
    const std::complex<double> val_b = mat_b.values[static_cast<size_t>(idx_b)];

    const int a_start = mat_a.col_index[static_cast<size_t>(row_idx_b)];
    const int a_end = mat_a.col_index[static_cast<size_t>(row_idx_b) + 1];

    for (int idx_a = a_start; idx_a < a_end; ++idx_a) {
      const int row_idx_a = mat_a.row_index[static_cast<size_t>(idx_a)];
      if (row_marker[static_cast<size_t>(row_idx_a)] != col_idx_b) {
        row_marker[static_cast<size_t>(row_idx_a)] = col_idx_b;
        active_rows.push_back(row_idx_a);
        accumulator[static_cast<size_t>(row_idx_a)] = mat_a.values[static_cast<size_t>(idx_a)] * val_b;
      } else {
        accumulator[static_cast<size_t>(row_idx_a)] += mat_a.values[static_cast<size_t>(idx_a)] * val_b;
      }
    }
  }

  std::ranges::sort(active_rows);

  for (const int row_idx_res : active_rows) {
    const std::complex<double> final_val = accumulator[static_cast<size_t>(row_idx_res)];
    if (IsValueNonZero(final_val)) {
      out_values.push_back(final_val);
      out_rows.push_back(row_idx_res);
    }
  }
}
}  // namespace

LiulinYComplexCcsOmp::LiulinYComplexCcsOmp(const InType &in) : BaseTask() {
  SetTypeOfTask(GetStaticTypeOfTask());
  GetInput() = in;
}

bool LiulinYComplexCcsOmp::ValidationImpl() {
  const auto &mat_a = GetInput().first;
  const auto &mat_b = GetInput().second;
  return mat_a.count_cols == mat_b.count_rows;
}

bool LiulinYComplexCcsOmp::PreProcessingImpl() {
  const auto &mat_a = GetInput().first;
  const auto &mat_b = GetInput().second;

  auto &mat_res = GetOutput();
  mat_res.count_rows = mat_a.count_rows;
  mat_res.count_cols = mat_b.count_cols;
  mat_res.values.clear();
  mat_res.row_index.clear();
  mat_res.col_index.assign(static_cast<size_t>(mat_res.count_cols) + 1, 0);

  return true;
}

bool LiulinYComplexCcsOmp::RunImpl() {
  const auto &mat_a = GetInput().first;
  const auto &mat_b = GetInput().second;
  auto &mat_res = GetOutput();

  const int rows_a_count = mat_a.count_rows;
  const int cols_b_count = mat_b.count_cols;

  std::vector<std::vector<std::complex<double>>> thread_values(static_cast<size_t>(cols_b_count));
  std::vector<std::vector<int>> thread_row_indices(static_cast<size_t>(cols_b_count));

#pragma omp parallel default(none) shared(mat_a, mat_b, thread_values, thread_row_indices, cols_b_count, rows_a_count)
  {
    std::vector<std::complex<double>> accumulator(static_cast<size_t>(rows_a_count), {0.0, 0.0});
    std::vector<int> active_rows;
    std::vector<int> row_marker(static_cast<size_t>(rows_a_count), -1);

#pragma omp for schedule(dynamic)
    for (int col_idx_b = 0; col_idx_b < cols_b_count; ++col_idx_b) {
      ProcessSingleColumn(col_idx_b, mat_a, mat_b, accumulator, active_rows, row_marker,
                          thread_values[static_cast<size_t>(col_idx_b)],
                          thread_row_indices[static_cast<size_t>(col_idx_b)]);
      active_rows.clear();
    }
  }

  mat_res.col_index[0] = 0;
  for (int col_idx = 0; col_idx < cols_b_count; ++col_idx) {
    const auto u_col_idx = static_cast<size_t>(col_idx);
    mat_res.values.insert(mat_res.values.end(), thread_values[u_col_idx].begin(), thread_values[u_col_idx].end());
    mat_res.row_index.insert(mat_res.row_index.end(), thread_row_indices[u_col_idx].begin(),
                             thread_row_indices[u_col_idx].end());
    mat_res.col_index[u_col_idx + 1] = static_cast<int>(mat_res.values.size());
  }

  return true;
}

bool LiulinYComplexCcsOmp::PostProcessingImpl() {
  return true;
}

}  // namespace liulin_y_complex_ccs
