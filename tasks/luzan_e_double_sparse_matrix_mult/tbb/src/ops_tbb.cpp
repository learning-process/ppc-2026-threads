#include "luzan_e_double_sparse_matrix_mult/tbb/include/ops_tbb.hpp"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <vector>

// #include <tbb/parallel_for.h>
// #include <tbb/blocked_range.h>
#include <tbb/tbb.h>

#include "luzan_e_double_sparse_matrix_mult/common/include/common.hpp"

namespace luzan_e_double_sparse_matrix_mult {
void LuzanEDoubleSparseMatrixMultTBB::MultiplyColumn(const SparseMatrix &a, const SparseMatrix &b, unsigned b_col,
                                                     std::vector<double> &tmp_col) {
  std::ranges::fill(tmp_col, 0.0);

  unsigned b_start = b.col_index[b_col];
  unsigned b_end = b.col_index[b_col + 1];

  for (unsigned b_pos = b_start; b_pos < b_end; ++b_pos) {
    double b_val = b.value[b_pos];
    unsigned k = b.row[b_pos];

    for (unsigned a_pos = a.col_index[k]; a_pos < a.col_index[k + 1]; ++a_pos) {
      tmp_col[a.row[a_pos]] += a.value[a_pos] * b_val;
    }
  }
}

void LuzanEDoubleSparseMatrixMultTBB::CompressColumn(const std::vector<double> &tmp_col, std::vector<double> &values,
                                                     std::vector<unsigned> &rows) {
  for (unsigned i = 0; i < tmp_col.size(); ++i) {
    if (fabs(tmp_col[i]) > kEPS) {
      values.push_back(tmp_col[i]);
      rows.push_back(i);
    }
  }
}

SparseMatrix LuzanEDoubleSparseMatrixMultTBB::CalcProdTBB(const SparseMatrix &a, const SparseMatrix &b) {
  SparseMatrix c(a.rows, b.cols);

  std::vector<std::vector<double>> values_per_col(b.cols);
  std::vector<std::vector<unsigned>> rows_per_col(b.cols);

  tbb::enumerable_thread_specific<std::vector<double>> tls_tmp([&] { return std::vector<double>(a.rows, 0.0); });

  tbb::parallel_for(tbb::blocked_range<unsigned>(0, b.cols), [&](const tbb::blocked_range<unsigned> &range) {
    for (unsigned b_col = range.begin(); b_col < range.end(); ++b_col) {
      auto &tmp_col = tls_tmp.local();

      MultiplyColumn(a, b, b_col, tmp_col);
      CompressColumn(tmp_col, values_per_col[b_col], rows_per_col[b_col]);
    }
  });

  BuildResult(c, values_per_col, rows_per_col);

  return c;
}

void LuzanEDoubleSparseMatrixMultTBB::BuildResult(SparseMatrix &c,
                                                  const std::vector<std::vector<double>> &values_per_col,
                                                  const std::vector<std::vector<unsigned>> &rows_per_col) {
  c.col_index.push_back(0);

  for (unsigned j = 0; j < values_per_col.size(); ++j) {
    for (size_t k = 0; k < values_per_col[j].size(); ++k) {
      c.value.push_back(values_per_col[j][k]);
      c.row.push_back(rows_per_col[j][k]);
    }
    c.col_index.push_back(c.value.size());
  }
}

LuzanEDoubleSparseMatrixMultTBB::LuzanEDoubleSparseMatrixMultTBB(const InType &in) {
  SetTypeOfTask(GetStaticTypeOfTask());
  GetInput() = in;
  // GetOutput() = 0;
}

bool LuzanEDoubleSparseMatrixMultTBB::ValidationImpl() {
  const auto &a = std::get<0>(GetInput());
  const auto &b = std::get<1>(GetInput());
  return a.GetCols() == b.GetRows() && a.GetCols() != 0 && a.GetRows() != 0 && b.GetCols() != 0;
}

bool LuzanEDoubleSparseMatrixMultTBB::PreProcessingImpl() {
  return true;
}

bool LuzanEDoubleSparseMatrixMultTBB::RunImpl() {
  const auto &a = std::get<0>(GetInput());
  const auto &b = std::get<1>(GetInput());

  GetOutput() = CalcProdTBB(a, b);
  return true;
}

bool LuzanEDoubleSparseMatrixMultTBB::PostProcessingImpl() {
  return true;
}

}  // namespace luzan_e_double_sparse_matrix_mult
