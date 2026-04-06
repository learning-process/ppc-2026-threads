#include "sizov_d_sparse_crs_mult/seq/include/ops_seq.hpp"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <vector>

#include "sizov_d_sparse_crs_mult/common/include/common.hpp"

namespace sizov_d_sparse_crs_mult {

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

SizovDSparseCRSMultSEQ::SizovDSparseCRSMultSEQ(const InType &in) {
  SetTypeOfTask(GetStaticTypeOfTask());
  GetInput() = in;
}

bool SizovDSparseCRSMultSEQ::ValidationImpl() {
  const auto &a = std::get<0>(GetInput());
  const auto &b = std::get<1>(GetInput());
  return IsValidCRS(a) && IsValidCRS(b) && a.cols == b.rows;
}

bool SizovDSparseCRSMultSEQ::PreProcessingImpl() {
  GetOutput() = {};
  return true;
}

bool SizovDSparseCRSMultSEQ::RunImpl() {
  const auto &a = std::get<0>(GetInput());
  const auto &b = std::get<1>(GetInput());
  auto &c = GetOutput();

  c.rows = a.rows;
  c.cols = b.cols;
  c.values.clear();
  c.col_indices.clear();
  c.row_ptr.assign(c.rows + 1, 0);

  std::vector<double> accum(c.cols, 0.0);
  std::vector<unsigned char> touched_flag(c.cols, 0);
  std::vector<std::size_t> touched_cols;

  for (std::size_t row_a = 0; row_a < a.rows; ++row_a) {
    for (std::size_t pos_a = a.row_ptr[row_a]; pos_a < a.row_ptr[row_a + 1]; ++pos_a) {
      const std::size_t row_b = a.col_indices[pos_a];
      const double val_a = a.values[pos_a];

      for (std::size_t pos_b = b.row_ptr[row_b]; pos_b < b.row_ptr[row_b + 1]; ++pos_b) {
        const std::size_t col_b = b.col_indices[pos_b];
        if (touched_flag[col_b] == 0U) {
          touched_flag[col_b] = 1U;
          touched_cols.push_back(col_b);
        }
        accum[col_b] += val_a * b.values[pos_b];
      }
    }

    std::ranges::sort(touched_cols);
    for (std::size_t col : touched_cols) {
      if (std::abs(accum[col]) > 1e-12) {
        c.col_indices.push_back(col);
        c.values.push_back(accum[col]);
      }
      accum[col] = 0.0;
      touched_flag[col] = 0U;
    }
    touched_cols.clear();
    c.row_ptr[row_a + 1] = c.values.size();
  }

  return true;
}

bool SizovDSparseCRSMultSEQ::PostProcessingImpl() {
  return true;
}

}  // namespace sizov_d_sparse_crs_mult
