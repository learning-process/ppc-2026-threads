#include "kapanova_s_sparse_matrix_mult_ccs_seq/seq/include/ops_seq.hpp"

#include <algorithm>
#include <cstddef>
#include <vector>

#include "kapanova_s_sparse_matrix_mult_ccs_seq/common/include/common.hpp"

namespace kapanova_s_sparse_matrix_mult_ccs_seq {

KapanovaSSparseMatrixMultCCSSeq::KapanovaSSparseMatrixMultCCSSeq(const InType &in) {
  SetTypeOfTask(GetStaticTypeOfTask());
  GetInput() = in;
}

bool KapanovaSSparseMatrixMultCCSSeq::ValidationImpl() {
  const auto &a = std::get<0>(GetInput());
  const auto &b = std::get<1>(GetInput());

  if (a.cols != b.rows) {
    return false;
  }
  if (a.rows == 0 || a.cols == 0 || b.rows == 0 || b.cols == 0) {
    return false;
  }
  if (a.col_ptrs.size() != a.cols + 1) {
    return false;
  }
  if (b.col_ptrs.size() != b.cols + 1) {
    return false;
  }
  if (a.col_ptrs[0] != 0 || b.col_ptrs[0] != 0) {
    return false;
  }
  if (a.col_ptrs[a.cols] != a.nnz) {
    return false;
  }
  if (b.col_ptrs[b.cols] != b.nnz) {
    return false;
  }
  if (a.values.size() != a.nnz || a.row_indices.size() != a.nnz) {
    return false;
  }
  if (b.values.size() != b.nnz || b.row_indices.size() != b.nnz) {
    return false;
  }

  return true;
}

bool KapanovaSSparseMatrixMultCCSSeq::PreProcessingImpl() {
  return true;
}

bool KapanovaSSparseMatrixMultCCSSeq::RunImpl() {
  const auto &a = std::get<0>(GetInput());
  const auto &b = std::get<1>(GetInput());
  OutType &c = GetOutput();

  c.values.clear();
  c.row_indices.clear();
  c.rows = a.rows;
  c.cols = b.cols;
  c.col_ptrs.assign(c.cols + 1, 0);

  std::vector<double> accum(a.rows, 0.0);
  std::vector<bool> nz_elem_rows(a.rows, false);
  std::vector<size_t> nnz_rows;

  for (size_t j = 0; j < b.cols; ++j) {
    c.col_ptrs[j] = c.values.size();

    for (size_t k = b.col_ptrs[j]; k < b.col_ptrs[j + 1]; ++k) {
      size_t ind = b.row_indices[k];
      double b_val = b.values[k];

      for (size_t zc = a.col_ptrs[ind]; zc < a.col_ptrs[ind + 1]; ++zc) {
        size_t i = a.row_indices[zc];
        double a_val = a.values[zc];

        accum[i] += a_val * b_val;
        if (!nz_elem_rows[i]) {
          nz_elem_rows[i] = true;
          nnz_rows.push_back(i);
        }
      }
    }

    std::ranges::sort(nnz_rows);

    for (size_t i : nnz_rows) {
      if (accum[i] != 0.0) {
        c.row_indices.push_back(i);
        c.values.push_back(accum[i]);
      }
      accum[i] = 0.0;
      nz_elem_rows[i] = false;
    }
    nnz_rows.clear();
  }

  c.nnz = c.values.size();
  c.col_ptrs[c.cols] = c.nnz;

  return true;
}

bool KapanovaSSparseMatrixMultCCSSeq::PostProcessingImpl() {
  return true;
}

}  // namespace kapanova_s_sparse_matrix_mult_ccs_seq
