#include "safaryan_a_sparse_matrix_mult_crs_seq/seq/include/ops_seq.hpp"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <vector>

#include "safaryan_a_sparse_matrix_mult_crs_seq/common/include/common.hpp"

namespace safaryan_a_sparse_matrix_mult_crs_seq {

SafaryanASparseMatrixMultCRSSeq::SafaryanASparseMatrixMultCRSSeq(const InType &in) {
  SetTypeOfTask(GetStaticTypeOfTask());
  GetInput() = in;
}

bool SafaryanASparseMatrixMultCRSSeq::IsMatrixValid(const CRSMatrix &m) {
  if (m.rows == 0 || m.cols == 0) {
    return false;
  }

  if (m.row_ptr.size() != m.rows + 1) {
    return false;
  }

  if (m.row_ptr.empty() || m.row_ptr[0] != 0) {
    return false;
  }

  if (m.values.size() != m.col_indices.size()) {
    return false;
  }
  if (m.nnz != m.values.size()) {
    return false;
  }

  if (m.row_ptr[m.rows] != m.nnz) {
    return false;
  }

  for (size_t i = 0; i < m.rows; ++i) {
    if (m.row_ptr[i] > m.row_ptr[i + 1]) {
      return false;
    }
    if (m.row_ptr[i + 1] > m.nnz) {
      return false;
    }
  }

  for (size_t idx = 0; idx < m.col_indices.size(); ++idx) {
    if (m.col_indices[idx] >= m.cols) {
      return false;
    }
  }

  return true;
}

bool SafaryanASparseMatrixMultCRSSeq::ValidationImpl() {
  const auto &a = std::get<0>(GetInput());
  const auto &b = std::get<1>(GetInput());

  if (!IsMatrixValid(a) || !IsMatrixValid(b)) {
    return false;
  }
  if (a.cols != b.rows) {
    return false;
  }

  return true;
}

bool SafaryanASparseMatrixMultCRSSeq::PreProcessingImpl() {
  const auto &a = std::get<0>(GetInput());
  const auto &b = std::get<1>(GetInput());

  OutType &c = GetOutput();
  c.values.clear();
  c.col_indices.clear();
  c.row_ptr.assign(a.rows + 1, 0);

  c.rows = a.rows;
  c.cols = b.cols;
  c.nnz = 0;

  return true;
}

bool SafaryanASparseMatrixMultCRSSeq::RunImpl() {
  const auto &a = std::get<0>(GetInput());
  const auto &b = std::get<1>(GetInput());
  OutType &c = GetOutput();

  c.values.clear();
  c.col_indices.clear();
  std::fill(c.row_ptr.begin(), c.row_ptr.end(), 0);

  std::vector<double> accum(b.cols, 0.0);
  std::vector<bool> used(b.cols, false);
  std::vector<size_t> used_cols;
  used_cols.reserve(64);

  constexpr double kEps = 1e-12;

  for (size_t i = 0; i < a.rows; ++i) {
    c.row_ptr[i] = c.values.size();

    for (size_t a_idx = a.row_ptr[i]; a_idx < a.row_ptr[i + 1]; ++a_idx) {
      const size_t k = a.col_indices[a_idx];
      const double a_val = a.values[a_idx];

      for (size_t b_idx = b.row_ptr[k]; b_idx < b.row_ptr[k + 1]; ++b_idx) {
        const size_t j = b.col_indices[b_idx];
        const double b_val = b.values[b_idx];

        accum[j] += a_val * b_val;

        if (!used[j]) {
          used[j] = true;
          used_cols.push_back(j);
        }
      }
    }

    std::sort(used_cols.begin(), used_cols.end());

    for (size_t j : used_cols) {
      const double v = accum[j];
      if (std::abs(v) > kEps) {
        c.col_indices.push_back(j);
        c.values.push_back(v);
      }
      accum[j] = 0.0;
      used[j] = false;
    }
    used_cols.clear();
  }

  c.nnz = c.values.size();
  c.row_ptr[c.rows] = c.nnz;

  return true;
}

bool SafaryanASparseMatrixMultCRSSeq::PostProcessingImpl() {
  return true;
}

}  // namespace safaryan_a_sparse_matrix_mult_crs_seq
