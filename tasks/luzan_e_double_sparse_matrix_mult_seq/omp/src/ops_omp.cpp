#include "luzan_e_double_sparse_matrix_mult_seq/omp/include/ops_omp.hpp"
// #include "luzan_e_double_sparse_matrix_mult_seq/common/include/common.hpp"
// #include "util/include/util.hpp"

#include <omp.h>

namespace luzan_e_double_sparse_matrix_mult_seq {

SparseMatrix SparseMatrix::CalcProdOMP(const SparseMatrix &b) const {
  SparseMatrix c(rows_, b.cols_);

  /// tmp storage
  std::vector<std::vector<double>> values_per_col(b.cols_);
  std::vector<std::vector<unsigned>> rows_per_col(b.cols_);

#pragma omp parallel for shared(b, values_per_col, rows_per_col) schedule(static) default(none) num_threads(2)
  for (int b_col = 0; b_col < static_cast<int>(b.cols_); b_col++) {
    std::vector<double> tmp_col(rows_, 0.0);

    unsigned b_rows_start = b.col_index_[b_col];
    unsigned b_rows_end = b.col_index_[b_col + 1];

    for (unsigned b_pos = b_rows_start; b_pos < b_rows_end; b_pos++) {
      double b_val = b.value_[b_pos];
      unsigned b_row = b.row_[b_pos];

      unsigned a_rows_start = col_index_[b_row];
      unsigned a_rows_end = col_index_[b_row + 1];

      for (unsigned a_pos = a_rows_start; a_pos < a_rows_end; a_pos++) {
        double a_val = value_[a_pos];
        unsigned a_row = row_[a_pos];

        tmp_col[a_row] += a_val * b_val;
      }
    }

    for (unsigned i = 0; i < rows_; i++) {
      if (fabs(tmp_col[i]) > kEPS) {
        values_per_col[b_col].push_back(tmp_col[i]);
        rows_per_col[b_col].push_back(i);
      }
    }
  }

  c.col_index_.push_back(0);
  for (unsigned j = 0; j < b.cols_; j++) {
    for (size_t k = 0; k < values_per_col[j].size(); k++) {
      c.value_.push_back(values_per_col[j][k]);
      c.row_.push_back(rows_per_col[j][k]);
    }
    c.col_index_.push_back(c.value_.size());
  }

  return c;
}

LuzanEDoubleSparseMatrixMultOMP::LuzanEDoubleSparseMatrixMultOMP(const InType &in) {
  SetTypeOfTask(GetStaticTypeOfTask());
  GetInput() = in;
  // GetOutput() = 0;
}

bool LuzanEDoubleSparseMatrixMultOMP::ValidationImpl() {
  const auto &a = std::get<0>(GetInput());
  const auto &b = std::get<1>(GetInput());
  return a.GetCols() == b.GetRows() && a.GetCols() != 0 && a.GetRows() != 0 && b.GetCols() != 0;
}

bool LuzanEDoubleSparseMatrixMultOMP::PreProcessingImpl() {
  return true;
}

bool LuzanEDoubleSparseMatrixMultOMP::RunImpl() {
  const auto &a = std::get<0>(GetInput());
  const auto &b = std::get<1>(GetInput());

  GetOutput() = a.CalcProdOMP(b);
  return true;
}

bool LuzanEDoubleSparseMatrixMultOMP::PostProcessingImpl() {
  return true;
}

}  // namespace luzan_e_double_sparse_matrix_mult_seq
