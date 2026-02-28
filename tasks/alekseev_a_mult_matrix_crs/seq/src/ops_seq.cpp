#include "alekseev_a_mult_matrix_crs/seq/include/ops_seq.hpp"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <vector>

namespace alekseev_a_mult_matrix_crs {

namespace {
bool IsValidCRS(const CRSMatrix &m) {
  if (m.rows == 0 || m.cols == 0) {
    return false;
  }
  if (m.row_ptr.size() != m.rows + 1) {
    return false;
  }
  if (m.row_ptr.empty() || m.row_ptr.front() != 0) {
    return false;
  }
  if (m.row_ptr.back() != m.values.size() || m.col_indices.size() != m.values.size()) {
    return false;
  }

  for (std::size_t i = 0; i < m.rows; ++i) {
    if (m.row_ptr[i] > m.row_ptr[i + 1]) {
      return false;
    }
  }
  return std::ranges::all_of(m.col_indices, [&m](std::size_t idx) { return idx < m.cols; });
}
}  // namespace

AlekseevAMultMatrixCRSSEQ::AlekseevAMultMatrixCRSSEQ(const InType &in) {
  SetTypeOfTask(GetStaticTypeOfTask());
  GetInput() = in;
}

bool AlekseevAMultMatrixCRSSEQ::ValidationImpl() {
  const auto &a = std::get<0>(GetInput());
  const auto &b = std::get<1>(GetInput());
  return IsValidCRS(a) && IsValidCRS(b) && a.cols == b.rows;
}

bool AlekseevAMultMatrixCRSSEQ::PreProcessingImpl() {
  GetOutput() = {};
  return true;
}

bool AlekseevAMultMatrixCRSSEQ::RunImpl() {
  const auto &a = std::get<0>(GetInput());
  const auto &b = std::get<1>(GetInput());
  auto &c = GetOutput();

  c.rows = a.rows;
  c.cols = b.cols;
  c.row_ptr.assign(c.rows + 1, 0);

  std::vector<double> accum(c.cols, 0.0);
  std::vector<int> touched_flag(c.cols, -1);
  std::vector<std::size_t> touched_cols;
  touched_cols.reserve(c.cols);

  for (std::size_t i = 0; i < a.rows; ++i) {
    for (std::size_t pos_a = a.row_ptr[i]; pos_a < a.row_ptr[i + 1]; ++pos_a) {
      std::size_t k = a.col_indices[pos_a];
      double val_a = a.values[pos_a];

      for (std::size_t pos_b = b.row_ptr[k]; pos_b < b.row_ptr[k + 1]; ++pos_b) {
        std::size_t j = b.col_indices[pos_b];

        if (touched_flag[j] != static_cast<int>(i)) {
          touched_flag[j] = static_cast<int>(i);
          touched_cols.push_back(j);
          accum[j] = 0.0;
        }
        accum[j] += val_a * b.values[pos_b];
      }
    }

    std::sort(touched_cols.begin(), touched_cols.end());

    for (auto col : touched_cols) {
      if (std::abs(accum[col]) > 1e-15) {
        c.values.push_back(accum[col]);
        c.col_indices.push_back(col);
      }
    }
    c.row_ptr[i + 1] = c.values.size();
    touched_cols.clear();
  }
  return true;
}

bool AlekseevAMultMatrixCRSSEQ::PostProcessingImpl() {
  return true;
}

}  // namespace alekseev_a_mult_matrix_crs
