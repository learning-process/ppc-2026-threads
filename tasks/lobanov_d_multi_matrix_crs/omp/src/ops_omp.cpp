#include "lobanov_d_multi_matrix_crs/omp/include/ops_omp.hpp"

#include <cmath>
#include <cstddef>
#include <map>
#include <vector>

#include "lobanov_d_multi_matrix_crs/common/include/common.hpp"
#include "util/include/util.hpp"

namespace lobanov_d_multi_matrix_crs {

LobanovMultyMatrixOMP::LobanovMultyMatrixOMP(const InType &in) {
  SetTypeOfTask(GetStaticTypeOfTask());
  GetInput() = in;
  GetOutput() = CompressedRowMatrix{};
}

bool LobanovMultyMatrixOMP::ValidationImpl() {
  const auto &matrix_a = GetInput().first;
  const auto &matrix_b = GetInput().second;

  if (matrix_a.row_count <= 0 || matrix_b.row_count <= 0 || matrix_a.column_count <= 0 || matrix_b.column_count <= 0) {
    return false;
  }
  if (matrix_a.column_count != matrix_b.row_count) {
    return false;
  }
  if (matrix_a.row_pointer_data.size() != static_cast<size_t>(matrix_a.row_count) + 1 ||
      matrix_b.row_pointer_data.size() != static_cast<size_t>(matrix_b.row_count) + 1) {
    return false;
  }
  if (static_cast<size_t>(matrix_a.non_zero_count) != matrix_a.value_data.size() ||
      static_cast<size_t>(matrix_a.non_zero_count) != matrix_a.column_index_data.size() ||
      static_cast<size_t>(matrix_b.non_zero_count) != matrix_b.value_data.size() ||
      static_cast<size_t>(matrix_b.non_zero_count) != matrix_b.column_index_data.size()) {
    return false;
  }
  return true;
}

bool LobanovMultyMatrixOMP::PreProcessingImpl() {
  const auto &matrix_a = GetInput().first;
  const auto &matrix_b = GetInput().second;

  auto &result = GetOutput();
  result.row_count = matrix_a.row_count;
  result.column_count = matrix_b.column_count;
  result.non_zero_count = 0;
  result.value_data.clear();
  result.column_index_data.clear();
  result.row_pointer_data.assign(static_cast<size_t>(result.row_count) + 1, 0);
  return true;
}

bool LobanovMultyMatrixOMP::RunImpl() {
  const auto &matrix_a = GetInput().first;
  const auto &matrix_b = GetInput().second;
  auto &result = GetOutput();

  const int rows_a = matrix_a.row_count;

  std::vector<std::map<int, double>> row_results(static_cast<size_t>(rows_a));

#pragma omp parallel for default(none) shared(matrix_a, matrix_b, row_results, rows_a) \
    num_threads(ppc::util::GetNumThreads()) schedule(dynamic)
  for (int i = 0; i < rows_a; ++i) {
    const int a_start = matrix_a.row_pointer_data[static_cast<size_t>(i)];
    const int a_end = matrix_a.row_pointer_data[static_cast<size_t>(i) + 1];

    for (int a_idx = a_start; a_idx < a_end; ++a_idx) {
      const int k = matrix_a.column_index_data[static_cast<size_t>(a_idx)];
      const double a_val = matrix_a.value_data[static_cast<size_t>(a_idx)];

      if (k >= matrix_b.row_count) {
        continue;
      }

      const int b_start = matrix_b.row_pointer_data[static_cast<size_t>(k)];
      const int b_end = matrix_b.row_pointer_data[static_cast<size_t>(k) + 1];

      for (int b_idx = b_start; b_idx < b_end; ++b_idx) {
        const int j = matrix_b.column_index_data[static_cast<size_t>(b_idx)];
        const double b_val = matrix_b.value_data[static_cast<size_t>(b_idx)];

        row_results[static_cast<size_t>(i)][j] += a_val * b_val;
      }
    }
  }

  int offset = 0;
  result.row_pointer_data[0] = 0;
  for (int i = 0; i < rows_a; ++i) {
    const auto &row = row_results[static_cast<size_t>(i)];
    for (const auto &[col, val] : row) {
      if (std::abs(val) > 1e-12) {
        result.column_index_data.push_back(col);
        result.value_data.push_back(val);
        ++offset;
      }
    }
    result.row_pointer_data[static_cast<size_t>(i) + 1] = offset;
  }
  result.non_zero_count = static_cast<int>(result.value_data.size());

  return true;
}

bool LobanovMultyMatrixOMP::PostProcessingImpl() {
  return true;
}

}  // namespace lobanov_d_multi_matrix_crs
