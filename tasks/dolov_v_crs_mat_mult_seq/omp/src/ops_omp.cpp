#include "dolov_v_crs_mat_mult_seq/omp/include/ops_omp.hpp"

#include <omp.h>

#include <cmath>
#include <cstddef>

namespace dolov_v_crs_mat_mult_seq {

DolovVCrsMatMultOmp::DolovVCrsMatMultOmp(const InType &in) {
  SetTypeOfTask(GetStaticTypeOfTask());
  GetInput() = in;
}

bool DolovVCrsMatMultOmp::ValidationImpl() {
  const auto &input_data = GetInput();
  if (input_data.size() != 2) {
    return false;
  }
  const auto &matrix_a = input_data[0];
  const auto &matrix_b = input_data[1];

  return matrix_a.num_cols == matrix_b.num_rows && matrix_a.num_rows > 0 && matrix_b.num_cols > 0 &&
         matrix_a.row_pointers.size() == static_cast<std::size_t>(matrix_a.num_rows + 1) &&
         matrix_b.row_pointers.size() == static_cast<std::size_t>(matrix_b.num_rows + 1);
}

bool DolovVCrsMatMultOmp::PreProcessingImpl() {
  auto &result = GetOutput();
  result.num_rows = GetInput()[0].num_rows;
  result.num_cols = GetInput()[1].num_cols;
  result.row_pointers.clear();
  result.values.clear();
  result.col_indices.clear();
  return true;
}

bool DolovVCrsMatMultOmp::RunImpl() {
  const auto &matrix_a = GetInput()[0];
  const auto &matrix_b = GetInput()[1];
  auto &result = GetOutput();

  const int m = matrix_a.num_rows;
  const int n = matrix_b.num_cols;

  std::vector<RowData> thread_rows(m);

#pragma omp parallel default(none) shared(matrix_a, matrix_b, thread_rows, m, n)
  {
    std::vector<double> row_accum(n, 0.0);
    std::vector<int> row_mask(n, -1);
    std::vector<int> active_cols;
    active_cols.reserve(n);

#pragma omp for schedule(dynamic, 16)
    for (int i = 0; i < m; ++i) {
      const int start_a = matrix_a.row_pointers[i];
      const int end_a = matrix_a.row_pointers[i + 1];

      for (int p_a = start_a; p_a < end_a; ++p_a) {
        int col_a = matrix_a.col_indices[p_a];
        double val_a = matrix_a.values[p_a];

        const int start_b = matrix_b.row_pointers[col_a];
        const int end_b = matrix_b.row_pointers[col_a + 1];

        for (int p_b = start_b; p_b < end_b; ++p_b) {
          int col_b = matrix_b.col_indices[p_b];
          double val_b = matrix_b.values[p_b];

          if (row_mask[col_b] != i) {
            row_mask[col_b] = i;
            active_cols.push_back(col_b);
          }
          row_accum[col_b] += val_a * val_b;
        }
      }

      for (std::size_t x = 1; x < active_cols.size(); ++x) {
        int key = active_cols[x];
        int y = static_cast<int>(x) - 1;
        while (y >= 0 && active_cols[y] > key) {
          active_cols[y + 1] = active_cols[y];
          y = y - 1;
        }
        active_cols[y + 1] = key;
      }

      for (int col : active_cols) {
        if (std::abs(row_accum[col]) > 1e-15) {
          thread_rows[i].cols.push_back(col);
          thread_rows[i].vals.push_back(row_accum[col]);
        }
        row_accum[col] = 0.0;
      }
      active_cols.clear();
    }
  }

  result.row_pointers.assign(static_cast<std::size_t>(m) + 1, 0);
  for (int i = 0; i < m; ++i) {
    result.row_pointers[i + 1] = result.row_pointers[i] + static_cast<int>(thread_rows[i].cols.size());
  }

  const int total_nnz = result.row_pointers[m];
  result.col_indices.resize(total_nnz);
  result.values.resize(total_nnz);

#pragma omp parallel for default(none) shared(result, thread_rows, m) schedule(static)
  for (int i = 0; i < m; ++i) {
    const int offset = result.row_pointers[i];
    const int nnz_in_row = result.row_pointers[i + 1] - offset;
    for (int j = 0; j < nnz_in_row; ++j) {
      result.col_indices[offset + j] = thread_rows[i].cols[j];
      result.values[offset + j] = thread_rows[i].vals[j];
    }
  }

  return true;
}

bool DolovVCrsMatMultOmp::PostProcessingImpl() {
  return true;
}

}  // namespace dolov_v_crs_mat_mult_seq
