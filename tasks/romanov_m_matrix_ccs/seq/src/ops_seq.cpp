#include "romanov_m_matrix_ccs/seq/include/ops_seq.hpp"

#include <algorithm>
#include <cmath>

namespace romanov_m_matrix_ccs {

RomanovMMatrixCCSSeq::RomanovMMatrixCCSSeq(const InType &in) {
  SetTypeOfTask(GetStaticTypeOfTask());
  GetInput() = in;
}

bool RomanovMMatrixCCSSeq::ValidationImpl() {
  const auto &left = GetInput().first;
  const auto &right = GetInput().second;

  if (left.cols_num != right.rows_num) {
    return false;
  }
  if (left.rows_num == 0 || left.cols_num == 0 || right.cols_num == 0) {
    return false;
  }
  if (left.col_ptrs.size() != left.cols_num + 1 || right.col_ptrs.size() != right.cols_num + 1) {
    return false;
  }

  return true;
}

bool RomanovMMatrixCCSSeq::PreProcessingImpl() {
  GetOutput().vals.clear();
  GetOutput().row_inds.clear();
  return true;
}

bool RomanovMMatrixCCSSeq::RunImpl() {
  const auto &A = GetInput().first;
  const auto &B = GetInput().second;
  auto &C = GetOutput();

  C.rows_num = A.rows_num;
  C.cols_num = B.cols_num;
  C.col_ptrs.assign(C.cols_num + 1, 0);

  std::vector<double> accumulator(A.rows_num, 0.0);
  std::vector<size_t> active_rows;
  std::vector<bool> row_mask(A.rows_num, false);

  for (size_t j = 0; j < B.cols_num; ++j) {
    C.col_ptrs[j] = C.vals.size();

    for (size_t kb = B.col_ptrs[j]; kb < B.col_ptrs[j + 1]; ++kb) {
      size_t k = B.row_inds[kb];
      double v_b = B.vals[kb];

      for (size_t ka = A.col_ptrs[k]; ka < A.col_ptrs[k + 1]; ++ka) {
        size_t i = A.row_inds[ka];
        if (!row_mask[i]) {
          row_mask[i] = true;
          active_rows.push_back(i);
        }
        accumulator[i] += A.vals[ka] * v_b;
      }
    }

    std::sort(active_rows.begin(), active_rows.end());

    for (size_t row_idx : active_rows) {
      if (std::abs(accumulator[row_idx]) > 1e-12) {
        C.vals.push_back(accumulator[row_idx]);
        C.row_inds.push_back(row_idx);
      }
      accumulator[row_idx] = 0.0;
      row_mask[row_idx] = false;
    }
    active_rows.clear();
  }

  C.nnz = C.vals.size();
  C.col_ptrs[C.cols_num] = C.nnz;
  return true;
}

bool RomanovMMatrixCCSSeq::PostProcessingImpl() {
  return true;
}

}  // namespace romanov_m_matrix_ccs
