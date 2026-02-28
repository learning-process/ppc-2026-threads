#include "shvetsova_k_mult_matrix_complex_col/seq/include/ops_seq.hpp"

#include <numeric>
#include <vector>

#include "shvetsova_k_mult_matrix_complex_col/common/include/common.hpp"
#include "util/include/util.hpp"

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
  const MatrixCCS &matrix_A = std::get<0>(GetInput());
  const MatrixCCS &matrix_B = std::get<1>(GetInput());
  auto &matrix_C = GetOutput();
  for (int i = 0; i < matrix_B.cols; i++) {
    std::vector<std::complex<double>> column_C(matrix_A.rows);
    for (int j = matrix_B.col_ptr[i]; j < matrix_B.col_ptr[i + 1]; j++) {
      int tmp_ind = matrix_B.row_ind[j];
      std::complex tmp_val = matrix_B.values[j];
      for (int t = matrix_A.col_ptr[tmp_ind]; t < matrix_A.col_ptr[tmp_ind + 1]; t++) {
        int row = matrix_A.row_ind[t];
        std::complex val_A = matrix_A.values[t];
        column_C[row] += tmp_val * val_A;
      }
    }
    int counter = 0;
    for (int k = 0; k < column_C.size(); k++) {
      if (column_C[k] != 0.0) {
        matrix_C.row_ind.push_back(k);
        matrix_C.values.push_back(column_C[k]);
        counter++;
      }
    }
    matrix_C.col_ptr.push_back(matrix_C.col_ptr[matrix_C.col_ptr.size() - 1] + counter);
  }
  return true;
}

bool ShvetsovaKMultMatrixComplexSEQ::PostProcessingImpl() {
  return true;
}

}  // namespace shvetsova_k_mult_matrix_complex_col
