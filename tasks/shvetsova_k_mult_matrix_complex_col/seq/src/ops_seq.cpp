#include "shvetsova_k_mult_matrix_complex_col/seq/include/ops_seq.hpp"

#include <complex>
#include <vector>

#include "shvetsova_k_mult_matrix_complex_col/common/include/common.hpp"

namespace shvetsova_k_mult_matrix_complex_col {

ShvetsovaKMultMatrixComplexSEQ::ShvetsovaKMultMatrixComplexSEQ(const InType &in) {
  SetTypeOfTask(GetStaticTypeOfTask());
  GetInput() = in;
  GetOutput() = MatrixCCS(0, 0, std::vector<int>{0}, std::vector<int>{}, std::vector<std::complex<double>>{});
}

bool ShvetsovaKMultMatrixComplexSEQ::ValidationImpl() {
  return ((std::get<0>(GetInput()).cols > 0) && (std::get<1>(GetInput()).cols > 0)) &&
         ((std::get<0>(GetInput()).rows > 0) && (std::get<1>(GetInput()).rows > 0));
}

bool ShvetsovaKMultMatrixComplexSEQ::PreProcessingImpl() {
  GetOutput().rows = std::get<0>(GetInput()).rows;
  GetOutput().cols = std::get<1>(GetInput()).cols;
  return GetOutput().rows == std::get<0>(GetInput()).rows && GetOutput().cols == std::get<1>(GetInput()).cols;
}

bool ShvetsovaKMultMatrixComplexSEQ::RunImpl() {
  const MatrixCCS &matrix_a = std::get<0>(GetInput());
  const MatrixCCS &matrix_b = std::get<1>(GetInput());
  auto &matrix_c = GetOutput();
  for (int i = 0; i < matrix_b.cols; i++) {
    std::vector<std::complex<double>> column_c(matrix_a.rows);
    for (int j = matrix_b.col_ptr[i]; j < matrix_b.col_ptr[i + 1]; j++) {
      int tmp_ind = matrix_b.row_ind[j];
      std::complex tmp_val = matrix_b.values[j];
      for (int ind = matrix_a.col_ptr[tmp_ind]; ind < matrix_a.col_ptr[tmp_ind + 1]; ind++) {
        int row = matrix_a.row_ind[ind];
        std::complex val_a = matrix_a.values[ind];
        column_c[row] += tmp_val * val_a;
      }
    }
    int counter = 0;
    for (size_t k = 0; k < column_c.size(); k++) {
      if (column_c[k] != 0.0) {
        matrix_c.row_ind.push_back(k);
        matrix_c.values.push_back(column_c[k]);
        counter++;
      }
    }
    matrix_c.col_ptr.push_back(matrix_c.col_ptr[static_cast<int>(matrix_c.col_ptr.size()) - 1] + counter);
  }
  return true;
}

bool ShvetsovaKMultMatrixComplexSEQ::PostProcessingImpl() {
  return true;
}

}  // namespace shvetsova_k_mult_matrix_complex_col
