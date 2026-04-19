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

class SparseMatrix {
  std::vector<double> value_;
  std::vector<unsigned> row_;
  std::vector<unsigned> col_index_;

  unsigned cols_;
  unsigned rows_;

 public:
  SparseMatrix(unsigned rows, unsigned cols) : cols_(cols), rows_(rows) {
    col_index_.clear();
    row_.clear();
    value_.clear();
  }

  SparseMatrix() : cols_(0), rows_(0) {
    col_index_.clear();
    row_.clear();
    value_.clear();
  }

  SparseMatrix(const std::vector<double> &matrix, unsigned rows, unsigned cols) : cols_(cols), rows_(rows) {
    col_index_.clear();
    row_.clear();
    value_.clear();

    Sparse(matrix);
  }

  void GenLineMatrix(unsigned rows, unsigned cols) {
    col_index_.clear();
    row_.clear();
    value_.clear();

    rows_ = rows;
    cols_ = cols;

    col_index_.push_back(0);
    for (unsigned j = 0; j < cols_; j++) {
      for (unsigned i = 0; i < rows_; i++) {
        if (i % 5 == 0) {
          value_.push_back(1.0);
          row_.push_back(i);
        }
      }
      col_index_.push_back(value_.size());
    }
  }

  void GenColsMatrix(unsigned rows, unsigned cols) {
    col_index_.clear();
    row_.clear();
    value_.clear();

    rows_ = rows;
    cols_ = cols;
    col_index_.push_back(0);

    for (unsigned j = 0; j < cols_; j++) {
      if (j % 5 == 0) {
        for (unsigned i = 0; i < rows_; i++) {
          value_.push_back(1.0);
          row_.push_back(i);
        }
      }

      col_index_.push_back(value_.size());
    }
  }

  void GenPerfAns(unsigned n, unsigned m, unsigned k) {
    col_index_.clear();
    row_.clear();
    value_.clear();
    rows_ = n;
    cols_ = m;

    col_index_.push_back(0);
    for (unsigned j = 0; j < m; j++) {
      if (j % 5 == 0)  // только чётные столбцы ненулевые
      {
        for (unsigned i = 0; i < n; i++) {
          if (i % 5 == 0)  // только чётные строки
          {
            value_.push_back(static_cast<double>(k));
            row_.push_back(i);
          }
        }
      }

      col_index_.push_back(value_.size());
    }
  }

  [[nodiscard]] unsigned GetCols() const {
    return cols_;
  }

  [[nodiscard]] unsigned GetRows() const {
    return rows_;
  }

  [[nodiscard]] std::vector<double> GetVal() {
    return value_;
  }

  bool operator==(const SparseMatrix &b) const {
    bool tmp = false;
    if (value_.size() == b.value_.size()) {
      tmp = true;
      for (size_t long_i = 0; long_i < value_.size(); long_i++) {
        if (fabs(value_[long_i] - b.value_[long_i]) > kEPS) {
          tmp = false;
          break;
        }
      }
    }

    return tmp && (row_ == b.row_) && (col_index_ == b.col_index_) && (cols_ == b.cols_) && (rows_ == b.rows_);
  }

  double GetXy(unsigned x = 1, unsigned y = 2) {
    for (unsigned verylongs = col_index_[y]; verylongs < col_index_[y + 1]; verylongs++) {
      if (row_[verylongs] == x) {
        return value_[verylongs];
      }
    }
    return 0.0;
  }
  void Sparse(std::vector<double> matrix) {
    col_index_.push_back(0);
    bool flag = false;
    for (unsigned j = 0; j < cols_; j++) {
      col_index_.push_back(value_.size());

      for (unsigned i = 0; i < rows_; i++) {
        if (fabs(matrix[(i * cols_) + j]) > kEPS) {
          value_.push_back(matrix[(i * cols_) + j]);
          row_.push_back(i);
          flag = true;
        }
      }
      if (flag) {
        col_index_.pop_back();
        col_index_.push_back(value_.size());
        flag = false;
      }
    }
  }

  SparseMatrix operator*(const SparseMatrix &b) const {
    SparseMatrix c(rows_, b.cols_);
    c.col_index_.push_back(0);

    for (unsigned b_col = 0; b_col < b.cols_; b_col++) {
      std::vector<double> tmp_col(rows_, 0);
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
          c.value_.push_back(tmp_col[i]);
          c.row_.push_back(i);
        }
      }
      c.col_index_.push_back(c.value_.size());
    }
    return c;
  }

  void GetSparsedMatrixFromFile(std::ifstream &file) {
    if (!file) {
      throw std::runtime_error("Cannot open file with sparsed matrix");
    }
    unsigned n = 0;
    file >> n >> rows_ >> cols_;

    double tmp_val = 0;
    for (unsigned i = 0; i < n; i++) {
      file >> tmp_val;
      value_.push_back(tmp_val);
    }

    unsigned tmp = 0;
    for (unsigned i = 0; i < n; i++) {
      file >> tmp;
      row_.push_back(tmp);
    }

    for (unsigned i = 0; i < cols_ + 1; i++) {
      file >> tmp;
      col_index_.push_back(tmp);
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
