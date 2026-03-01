#pragma once

#include <cmath>
#include <std::vector>
#include <string>
#include <tuple>

#include "task/include/task.hpp"

namespace luzan_e_double_sparse_matrix_mult_seq {

using InType = std::tuple<Sparse_matrix, Sparse_matrix>;
using OutType = Sparse_matrix;
using TestType = std::tuple<int, std::string>;
using BaseTask = ppc::task::Task<InType, OutType>;

const double EPS = 1e-8;

class Sparse_matrix {
  std::vector<double> value;
  std::vector<unsigned> row;
  std::vector<unsigned> col_index;

  unsigned cols;
  unsigned rows;

 public:
  Sparse_matrix(unsigned rows_, unsigned cols_) : cols(cols_), rows(rows_) {
    col_index.clear();
    row.clear();
    value.clear();
  }

  Sparse_matrix() : cols(0), rows(0) {
    col_index.clear();
    row.clear();
    value.clear();
  }

  Sparse_matrix(const std::vector<double> matrix, unsigned rows_, unsigned cols_) : cols(cols_), rows(rows_) {
    col_index.clear();
    row.clear();
    value.clear();

    sparse(matrix);
  }

  double getXY(int x = 1, int y = 2) {
    for (int s = col_index[y]; s < col_index[y + 1]; s++) {
      if (row[s] == x) {
        return value[s];
      }
    }
    return 0.0;
  }
  void sparse(std::vector<double> matrix) {
    col_index.push_back(0);
    bool flag = false;
    for (unsigned j = 0; j < cols; j++) {
      col_index.push_back(value.size());

      for (unsigned i = 0; i < rows; i++) {
        if (fabs(matrix[(i * cols) + j]) > EPS) {
          value.push_back(matrix[(i * cols) + j]);
          row.push_back(i);
          flag = true;
        }
      }
      if (flag == true)  // можно сравнивать последний из col_ind и value.size, но зачем
      {
        col_index.pop_back();
        col_index.push_back(value.size());
        flag = false;
      }
    }
  }

  Sparse_matrix operator*(const Sparse_matrix &B) const {
    Sparse_matrix C(rows, B.cols);
    C.col_index.push_back(0);

    for (unsigned b_col = 0; b_col < B.cols; b_col++) {
      std::vector<double> tmp_col(rows, 0);
      unsigned b_rows_start = B.col_index[b_col];
      unsigned b_rows_end = B.col_index[b_col + 1];

      for (unsigned b_pos = b_rows_start; b_pos < b_rows_end; b_pos++) {
        double b_val = B.value[b_pos];
        unsigned b_row = B.row[b_pos];

        unsigned a_rows_start = col_index[b_row];
        unsigned a_rows_end = col_index[b_row + 1];

        for (unsigned a_pos = a_rows_start; a_pos < a_rows_end; a_pos++) {
          double a_val = value[a_pos];
          unsigned a_row = row[a_pos];
          tmp_col[a_row] += a_val * b_val;
        }
      }
      for (unsigned i = 0; i < rows; i++) {
        if (fabs(tmp_col[i]) > EPS) {
          C.value.push_back(tmp_col[i]);
          C.row.push_back(i);
        }
      }
      C.col_index.push_back(C.value.size());
    }
    return C;
  }
};

}  // namespace luzan_e_double_sparse_matrix_mult_seq
