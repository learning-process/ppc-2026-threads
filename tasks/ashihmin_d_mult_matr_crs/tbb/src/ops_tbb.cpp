#include "ashihmin_d_mult_matr_crs/tbb/include/ops_tbb.hpp"

#include <tbb/tbb.h>

#include <cmath>
#include <map>
#include <vector>

#include "ashihmin_d_mult_matr_crs/common/include/common.hpp"

namespace ashihmin_d_mult_matr_crs {

AshihminDMultMatrCrsTBB::AshihminDMultMatrCrsTBB(const InType &in) {
  SetTypeOfTask(GetStaticTypeOfTask());
  GetInput() = in;
}

bool AshihminDMultMatrCrsTBB::ValidationImpl() {
  return GetInput().first.cols == GetInput().second.rows;
}

bool AshihminDMultMatrCrsTBB::PreProcessingImpl() {
  auto &matrix_c = GetOutput();
  matrix_c.rows = GetInput().first.rows;
  matrix_c.cols = GetInput().second.cols;
  matrix_c.row_ptr.assign(matrix_c.rows + 1, 0);
  matrix_c.values.clear();
  matrix_c.col_index.clear();
  return true;
}

bool AshihminDMultMatrCrsTBB::RunImpl() {
  const auto &matrix_a = GetInput().first;
  const auto &matrix_b = GetInput().second;
  auto &matrix_c = GetOutput();
  int rows_a = matrix_a.rows;

  std::vector<std::vector<int>> local_cols(rows_a);
  std::vector<std::vector<double>> local_vals(rows_a);

  tbb::parallel_for(0, rows_a, [&](int i) {
    std::map<int, double> row_accumulator;
    for (int j = matrix_a.row_ptr[i]; j < matrix_a.row_ptr[i + 1]; ++j) {
      int col_a = matrix_a.col_index[j];
      double val_a = matrix_a.values[j];
      for (int k = matrix_b.row_ptr[col_a]; k < matrix_b.row_ptr[col_a + 1]; ++k) {
        row_accumulator[matrix_b.col_index[k]] += val_a * matrix_b.values[k];
      }
    }
    for (const auto &[col, val] : row_accumulator) {
      if (std::abs(val) > 1e-15) {
        local_cols[i].push_back(col);
        local_vals[i].push_back(val);
      }
    }
  });

  for (int i = 0; i < rows_a; ++i) {
    matrix_c.col_index.insert(matrix_c.col_index.end(), local_cols[i].begin(), local_cols[i].end());
    matrix_c.values.insert(matrix_c.values.end(), local_vals[i].begin(), local_vals[i].end());
    matrix_c.row_ptr[i + 1] = static_cast<int>(matrix_c.values.size());
  }
  return true;
}

bool AshihminDMultMatrCrsTBB::PostProcessingImpl() {
  return true;
}

}  // namespace ashihmin_d_mult_matr_crs
