#pragma once

#include <complex>
#include <string>
#include <tuple>
#include <vector>

#include "task/include/task.hpp"

namespace shvetsova_k_mult_matrix_complex_col {

struct MatrixCCS {
  int rows = 0;
  int cols = 0;
  std::vector<int> col_ptr;
  std::vector<int> row_ind;
  std::vector<std::complex<double>> values;

  bool operator==(const MatrixCCS &matrix) const {
    return rows == matrix.rows && cols == matrix.cols && col_ptr == matrix.col_ptr && row_ind == matrix.row_ind &&
           values == matrix.values;
  }
};
using InType = std::tuple<MatrixCCS, MatrixCCS>;
using OutType = MatrixCCS;
using TestType = std::tuple<int, std::string>;
using BaseTask = ppc::task::Task<InType, OutType>;

}  // namespace shvetsova_k_mult_matrix_complex_col
