#pragma once

#include <cmath>
#include <cstddef>
#include <fstream>
#include <stdexcept>
#include <string>
#include <tuple>
#include <vector>

#include "task/include/task.hpp"

namespace luzan_e_double_sparse_matrix_mult {

const double kEPS = 1e-8;

struct SparseMatrix {
  std::vector<double> value;
  std::vector<unsigned> row;
  std::vector<unsigned> col_index;

  unsigned cols;
  unsigned rows;

  SparseMatrix(unsigned rows_out, unsigned cols_out) : cols(cols_out), rows(rows_out) {
    col_index.clear();
    row.clear();
    value.clear();
  }

  SparseMatrix() : cols(0), rows(0) {
    col_index.clear();
    row.clear();
    value.clear();
  }

  SparseMatrix(const std::vector<double> &matrix, unsigned rows_out, unsigned cols_out)
      : cols(cols_out), rows(rows_out) {
    col_index.clear();
    row.clear();
    value.clear();

    Sparse(matrix);
  }

  void GenLineMatrix(unsigned rows_out, unsigned cols_out) {
    col_index.clear();
    row.clear();
    value.clear();

    rows = rows_out;
    cols = cols_out;

    col_index.push_back(0);
    for (unsigned j = 0; j < cols; j++) {
      for (unsigned i = 0; i < rows; i++) {
        if (i % 5 == 0) {
          value.push_back(1.0);
          row.push_back(i);
        }
      }
      col_index.push_back(value.size());
    }
  }

  void GenColsMatrix(unsigned rows_out, unsigned cols_out) {
    col_index.clear();
    row.clear();
    value.clear();

    rows = rows_out;
    cols = cols_out;
    col_index.push_back(0);

    for (unsigned j = 0; j < cols; j++) {
      if (j % 5 == 0) {
        for (unsigned i = 0; i < rows; i++) {
          value.push_back(1.0);
          row.push_back(i);
        }
      }

      col_index.push_back(value.size());
    }
  }

  void GenPerfAns(unsigned n, unsigned m, unsigned k) {
    col_index.clear();
    row.clear();
    value.clear();
    rows = n;
    cols = m;

    col_index.push_back(0);
    for (unsigned j = 0; j < m; j++) {
      if (j % 5 == 0)  // только чётные столбцы ненулевые
      {
        for (unsigned i = 0; i < n; i++) {
          if (i % 5 == 0)  // только чётные строки
          {
            value.push_back(static_cast<double>(k));
            row.push_back(i);
          }
        }
      }

      col_index.push_back(value.size());
    }
  }

  [[nodiscard]] unsigned GetCols() const {
    return cols;
  }

  [[nodiscard]] unsigned GetRows() const {
    return rows;
  }

  [[nodiscard]] std::vector<double> GetVal() const {
    return value;
  }

  bool operator==(const SparseMatrix &b) const {
    bool tmp = false;
    if (value.size() == b.value.size()) {
      tmp = true;
      for (size_t long_i = 0; long_i < value.size(); long_i++) {
        if (fabs(value[long_i] - b.value[long_i]) > kEPS) {
          tmp = false;
          break;
        }
      }
    }

    return tmp && (row == b.row) && (col_index == b.col_index) && (cols == b.cols) && (rows == b.rows);
  }

  double GetXy(unsigned x = 1, unsigned y = 2) {
    for (unsigned verylongs = col_index[y]; verylongs < col_index[y + 1]; verylongs++) {
      if (row[verylongs] == x) {
        return value[verylongs];
      }
    }
    return 0.0;
  }
  void Sparse(const std::vector<double> &matrix) {
    col_index.push_back(0);
    bool flag = false;
    for (unsigned j = 0; j < cols; j++) {
      col_index.push_back(value.size());

      for (unsigned i = 0; i < rows; i++) {
        if (fabs(matrix[(i * cols) + j]) > kEPS) {
          value.push_back(matrix[(i * cols) + j]);
          row.push_back(i);
          flag = true;
        }
      }
      if (flag) {
        col_index.pop_back();
        col_index.push_back(value.size());
        flag = false;
      }
    }
  }

  SparseMatrix operator*(const SparseMatrix &b) const {
    SparseMatrix c(rows, b.cols);
    c.col_index.push_back(0);

    for (unsigned b_col = 0; b_col < b.cols; b_col++) {
      std::vector<double> tmp_col(rows, 0);
      unsigned b_rows_start = b.col_index[b_col];
      unsigned b_rows_end = b.col_index[b_col + 1];

      for (unsigned b_pos = b_rows_start; b_pos < b_rows_end; b_pos++) {
        double b_val = b.value[b_pos];
        unsigned b_row = b.row[b_pos];

        unsigned a_rows_start = col_index[b_row];
        unsigned a_rows_end = col_index[b_row + 1];

        for (unsigned a_pos = a_rows_start; a_pos < a_rows_end; a_pos++) {
          double a_val = value[a_pos];
          unsigned a_row = row[a_pos];
          tmp_col[a_row] += a_val * b_val;
        }
      }
      for (unsigned i = 0; i < rows; i++) {
        if (fabs(tmp_col[i]) > kEPS) {
          c.value.push_back(tmp_col[i]);
          c.row.push_back(i);
        }
      }
      c.col_index.push_back(c.value.size());
    }
    return c;
  }

  void GetSparsedMatrixFromFile(std::ifstream &file) {
    if (!file) {
      throw std::runtime_error("Cannot open file with sparsed matrix");
    }
    value.clear();
    row.clear();
    col_index.clear();
    unsigned n = 0;
    file >> n >> rows >> cols;

    double tmp_val = 0;
    for (unsigned i = 0; i < n; i++) {
      file >> tmp_val;
      value.push_back(tmp_val);
    }

    unsigned tmp = 0;
    for (unsigned i = 0; i < n; i++) {
      file >> tmp;
      row.push_back(tmp);
    }

    for (unsigned i = 0; i < cols + 1; i++) {
      file >> tmp;
      col_index.push_back(tmp);
    }
  }
};

inline SparseMatrix GetFromFile(std::ifstream &file) {
  size_t r = 0;
  size_t c = 0;
  file >> r >> c;

  std::vector<double> dense(r * c);

  for (unsigned i = 0; i < r; i++) {
    for (unsigned j = 0; j < c; j++) {
      file >> dense[(i * c) + j];
    }
  }
  SparseMatrix a(dense, r, c);
  return a;
};

using InType = std::tuple<SparseMatrix, SparseMatrix>;
using OutType = SparseMatrix;
using TestType = std::tuple<std::string, std::string>;
using BaseTask = ppc::task::Task<InType, OutType>;

}  // namespace luzan_e_double_sparse_matrix_mult
