#include "luzan_e_double_sparse_matrix_mult/omp/include/ops_omp.hpp"

#include <cmath>
#include <cstddef>
#include <vector>

#include "luzan_e_double_sparse_matrix_mult/common/include/common.hpp"

namespace luzan_e_double_sparse_matrix_mult {
SparseMatrix LuzanEDoubleSparseMatrixMultOMP::CalcProdOMP(const SparseMatrix &a, const SparseMatrix &b) {
  SparseMatrix c(a.rows, b.cols);

  /// tmp storage
  std::vector<std::vector<double>> values_per_col(b.cols);
  std::vector<std::vector<unsigned>> rowsper_col(b.cols);

#pragma omp parallel for shared(a, b, values_per_col, rowsper_col, kEPS) schedule(static) default(none)
  for (unsigned b_col = 0; b_col < static_cast<unsigned>(b.cols); b_col++) {
    std::vector<double> tmp_col(a.rows, 0.0);

    unsigned b_rowsstart = b.col_index[b_col];
    unsigned b_rowsend = b.col_index[b_col + 1];

    for (unsigned b_pos = b_rowsstart; b_pos < b_rowsend; b_pos++) {
      double b_val = b.value[b_pos];
      unsigned b_row = b.row[b_pos];

      unsigned a_rowsstart = a.col_index[b_row];
      unsigned a_rowsend = a.col_index[b_row + 1];

      for (unsigned a_pos = a_rowsstart; a_pos < a_rowsend; a_pos++) {
        double a_val = a.value[a_pos];
        unsigned a_row = a.row[a_pos];

        tmp_col[a_row] += a_val * b_val;
      }
    }

    for (unsigned i = 0; i < a.rows; i++) {
      if (fabs(tmp_col[i]) > kEPS) {
        values_per_col[b_col].push_back(tmp_col[i]);
        rowsper_col[b_col].push_back(i);
      }
    }
  }

  c.col_index.push_back(0);
  for (unsigned j = 0; j < b.cols; j++) {
    for (size_t k = 0; k < values_per_col[j].size(); k++) {
      c.value.push_back(values_per_col[j][k]);
      c.row.push_back(rowsper_col[j][k]);
    }
    c.col_index.push_back(c.value.size());
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

  GetOutput() = CalcProdOMP(a, b);
  return true;
}

bool LuzanEDoubleSparseMatrixMultOMP::PostProcessingImpl() {
  return true;
}

}  // namespace luzan_e_double_sparse_matrix_mult
