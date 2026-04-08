#include "sizov_d_sparse_crs_mult/omp/include/ops_omp.hpp"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <utility>
#include <vector>

#include "sizov_d_sparse_crs_mult/common/include/common.hpp"

namespace sizov_d_sparse_crs_mult {

namespace {

void AccumulateRowProducts(std::size_t row_idx, const CRSMatrix &a, const CRSMatrix &b, std::vector<double> &accum,
                           std::vector<unsigned char> &touched, std::vector<std::size_t> &touched_cols) {
  for (std::size_t a_idx = a.row_ptr[row_idx]; a_idx < a.row_ptr[row_idx + 1]; ++a_idx) {
    const std::size_t k = a.col_indices[a_idx];
    const double a_val = a.values[a_idx];

    for (std::size_t b_idx = b.row_ptr[k]; b_idx < b.row_ptr[k + 1]; ++b_idx) {
      const std::size_t j = b.col_indices[b_idx];
      if (touched[j] == 0U) {
        touched[j] = 1U;
        touched_cols.push_back(j);
      }
      accum[j] += a_val * b.values[b_idx];
    }
  }
}

void FlushRowToEntries(std::vector<double> &accum, std::vector<unsigned char> &touched,
                       std::vector<std::size_t> &touched_cols, std::vector<std::pair<std::size_t, double>> &row) {
  std::ranges::sort(touched_cols);
  row.clear();
  row.reserve(touched_cols.size());
  for (std::size_t col : touched_cols) {
    const double value = accum[col];
    if (std::abs(value) > 1e-12) {
      row.emplace_back(col, value);
    }
    accum[col] = 0.0;
    touched[col] = 0U;
  }
  touched_cols.clear();
}

}  // namespace

SizovDSparseCRSMultOMP::SizovDSparseCRSMultOMP(const InType &in) {
  SetTypeOfTask(GetStaticTypeOfTask());
  GetInput() = in;
}

bool SizovDSparseCRSMultOMP::ValidationImpl() {
  const auto &input = GetInput();
  const auto &a = std::get<0>(input);
  const auto &b = std::get<1>(input);

  if (a.cols != b.rows) {
    return false;
  }
  if (a.row_ptr.size() != a.rows + 1) {
    return false;
  }
  if (b.row_ptr.size() != b.rows + 1) {
    return false;
  }
  if (a.values.size() != a.col_indices.size()) {
    return false;
  }
  if (b.values.size() != b.col_indices.size()) {
    return false;
  }
  if (a.row_ptr.back() != a.values.size()) {
    return false;
  }
  if (b.row_ptr.back() != b.values.size()) {
    return false;
  }

  return true;
}

bool SizovDSparseCRSMultOMP::PreProcessingImpl() {
  GetOutput() = CRSMatrix{};
  return true;
}

bool SizovDSparseCRSMultOMP::RunImpl() {
  const auto &input = GetInput();
  const auto &a = std::get<0>(input);
  const auto &b = std::get<1>(input);

  CRSMatrix c;
  c.rows = a.rows;
  c.cols = b.cols;

  std::vector<std::vector<std::pair<std::size_t, double>>> row_entries(c.rows);

#pragma omp parallel default(none) shared(a, b, c, row_entries)
  {
    std::vector<double> accum(c.cols, 0.0);
    std::vector<unsigned char> touched(c.cols, 0);
    std::vector<std::size_t> touched_cols;
    touched_cols.reserve(256);

    const auto rows_end = static_cast<int64_t>(a.rows);
#pragma omp for schedule(dynamic, 16)
    for (int64_t i = 0; i < rows_end; ++i) {
      const auto row_idx = static_cast<std::size_t>(i);
      auto &row = row_entries[row_idx];
      AccumulateRowProducts(row_idx, a, b, accum, touched, touched_cols);
      FlushRowToEntries(accum, touched, touched_cols, row);
    }
  }

  c.row_ptr.resize(c.rows + 1, 0);
  for (std::size_t i = 0; i < c.rows; ++i) {
    c.row_ptr[i + 1] = c.row_ptr[i] + row_entries[i].size();
  }

  c.values.resize(c.row_ptr.back());
  c.col_indices.resize(c.row_ptr.back());

#pragma omp parallel for default(none) shared(c, row_entries)
  for (std::size_t i = 0; i < c.rows; ++i) {
    std::size_t out_pos = c.row_ptr[i];
    for (const auto &[col, value] : row_entries[i]) {
      c.col_indices[out_pos] = col;
      c.values[out_pos] = value;
      ++out_pos;
    }
  }

  GetOutput() = std::move(c);
  return true;
}

bool SizovDSparseCRSMultOMP::PostProcessingImpl() {
  return true;
}

}  // namespace sizov_d_sparse_crs_mult
