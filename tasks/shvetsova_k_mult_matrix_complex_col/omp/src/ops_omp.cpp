#include "../../omp/include/ops_omp.hpp"

#include <complex>
#include <utility>
#include <vector>

#include "omp.h"
#include "shvetsova_k_mult_matrix_complex_col/common/include/common.hpp"

namespace shvetsova_k_mult_matrix_complex_col {

ShvetsovaKMultMatrixComplexOMP::ShvetsovaKMultMatrixComplexOMP(const InType &in) {
  SetTypeOfTask(GetStaticTypeOfTask());
  GetInput() = in;
  GetOutput() = MatrixCCS(0, 0, std::vector<int>{0}, std::vector<int>{}, std::vector<std::complex<double>>{});
}

bool ShvetsovaKMultMatrixComplexOMP::ValidationImpl() {
  return true;
}

bool ShvetsovaKMultMatrixComplexOMP::PreProcessingImpl() {
  const auto &matrix_a = std::get<0>(GetInput());
  const auto &matrix_b = std::get<1>(GetInput());

  auto &matrix_c = GetOutput();
  matrix_c.rows = matrix_a.rows;
  matrix_c.cols = matrix_b.cols;
  matrix_c.row_ind.clear();
  matrix_c.values.clear();
  matrix_c.col_ptr.clear();
  matrix_c.col_ptr.push_back(0);
  return true;
}

bool ShvetsovaKMultMatrixComplexOMP::RunImpl() {
  const MatrixCCS &matrix_a = std::get<0>(GetInput());
  const MatrixCCS &matrix_b = std::get<1>(GetInput());
  struct SparseColumn {
    std::vector<int> rows;
    std::vector<std::complex<double>> vals;
  };

  std::vector<SparseColumn> columns_c(matrix_b.cols);
  auto &matrix_c = GetOutput();

#pragma omp parallel default(none) shared(matrix_a, matrix_b, columns_c)
  {
    std::vector<std::complex<double>> column_c(matrix_a.rows, {0.0, 0.0});
#pragma omp for
    for (int i = 0; i < matrix_b.cols; i++) {
      for (auto &val : column_c) {
        val = {0.0, 0.0};
      }
      for (int j = matrix_b.col_ptr[i]; j < matrix_b.col_ptr[i + 1]; j++) {
        int tmp_ind = matrix_b.row_ind[j];
        std::complex tmp_val = matrix_b.values[j];
        for (int ind = matrix_a.col_ptr[tmp_ind]; ind < matrix_a.col_ptr[tmp_ind + 1]; ind++) {
          int row = matrix_a.row_ind[ind];
          std::complex val_a = matrix_a.values[ind];
          column_c[row] += tmp_val * val_a;
        }
      }

      for (int index = 0; std::cmp_less(index, column_c.size()); ++index) {
        if (column_c[index].real() != 0.0 || column_c[index].imag() != 0.0) {
          columns_c[i].rows.push_back(index);
          columns_c[i].vals.push_back(column_c[index]);
        }
      }
    }
  }

  for (int i = 0; i < matrix_b.cols; i++) {
    matrix_c.row_ind.insert(matrix_c.row_ind.end(), columns_c[i].rows.begin(), columns_c[i].rows.end());
    matrix_c.values.insert(matrix_c.values.end(), columns_c[i].vals.begin(), columns_c[i].vals.end());
    matrix_c.col_ptr.push_back(static_cast<int>(matrix_c.row_ind.size()));
  }

  return true;
}

bool ShvetsovaKMultMatrixComplexOMP::PostProcessingImpl() {
  return true;
}

}  // namespace shvetsova_k_mult_matrix_complex_col
