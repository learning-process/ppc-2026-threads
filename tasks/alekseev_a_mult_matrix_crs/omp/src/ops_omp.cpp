#include "alekseev_a_mult_matrix_crs/omp/include/ops_omp.hpp"

#include <omp.h>

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <utility>
#include <vector>

#include "alekseev_a_mult_matrix_crs/common/include/common.hpp"

namespace alekseev_a_mult_matrix_crs {

AlekseevAMultMatrixCRSOMP::AlekseevAMultMatrixCRSOMP(const InType &in) {
  SetTypeOfTask(GetStaticTypeOfTask());
  GetInput() = in;
}

bool AlekseevAMultMatrixCRSOMP::ValidationImpl() {
  const auto &a = std::get<0>(GetInput());
  const auto &b = std::get<1>(GetInput());
  return a.cols == b.rows && !a.row_ptr.empty() && !b.row_ptr.empty();
}

bool AlekseevAMultMatrixCRSOMP::PreProcessingImpl() {
  GetOutput() = {};
  return true;
}

void AlekseevAMultMatrixCRSOMP::ProcessRow(std::size_t i, const CRSMatrix &a, const CRSMatrix &b,
                                           std::vector<double> &temp_v, std::vector<std::size_t> &temp_c,
                                           std::vector<double> &accum, std::vector<int> &touched_flag,
                                           std::vector<std::size_t> &touched_cols) {
  for (std::size_t pos_a = a.row_ptr[i]; pos_a < a.row_ptr[i + 1]; ++pos_a) {
    std::size_t k = a.col_indices[pos_a];
    double val_a = a.values[pos_a];

    for (std::size_t pos_b = b.row_ptr[k]; pos_b < b.row_ptr[k + 1]; ++pos_b) {
      std::size_t j = b.col_indices[pos_b];

      if (std::cmp_not_equal(touched_flag[j], i)) {
        touched_flag[j] = static_cast<int>(i);
        touched_cols.push_back(j);
        accum[j] = 0.0;
      }
      accum[j] += val_a * b.values[pos_b];
    }
  }

  std::ranges::sort(touched_cols);

  for (auto col : touched_cols) {
    if (std::abs(accum[col]) > 1e-15) {
      temp_v.push_back(accum[col]);
      temp_c.push_back(col);
    }
  }
  touched_cols.clear();
}

bool AlekseevAMultMatrixCRSOMP::RunImpl() {
  const auto &a = std::get<0>(GetInput());
  const auto &b = std::get<1>(GetInput());
  auto &c = GetOutput();

  c.rows = a.rows;
  c.cols = b.cols;
  c.row_ptr.assign(c.rows + 1, 0);

  std::vector<std::vector<double>> temp_values(c.rows);
  std::vector<std::vector<std::size_t>> temp_cols(c.rows);

  const int rows_limit = static_cast<int>(a.rows);

#pragma omp parallel default(none) shared(a, b, c, temp_values, temp_cols, rows_limit)
  {
    std::vector<double> accum(c.cols, 0.0);
    std::vector<int> touched_flag(c.cols, -1);
    std::vector<std::size_t> touched_cols;
    touched_cols.reserve(c.cols);

#pragma omp for schedule(dynamic)
    for (int i = 0; i < rows_limit; ++i) {
      ProcessRow(static_cast<std::size_t>(i), a, b, temp_values[i], temp_cols[i], accum, touched_flag, touched_cols);
    }
  }

  for (std::size_t i = 0; i < c.rows; ++i) {
    c.row_ptr[i + 1] = c.row_ptr[i] + temp_values[i].size();
  }

  c.values.reserve(c.row_ptr[c.rows]);
  c.col_indices.reserve(c.row_ptr[c.rows]);

  for (std::size_t i = 0; i < c.rows; ++i) {
    c.values.insert(c.values.end(), temp_values[i].begin(), temp_values[i].end());
    c.col_indices.insert(c.col_indices.end(), temp_cols[i].begin(), temp_cols[i].end());
  }

  return true;
}

bool AlekseevAMultMatrixCRSOMP::PostProcessingImpl() {
  return true;
}

}  // namespace alekseev_a_mult_matrix_crs
