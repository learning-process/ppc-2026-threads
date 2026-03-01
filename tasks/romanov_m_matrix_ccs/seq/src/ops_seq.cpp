#include "romanov_m_matrix_ccs/seq/include/ops_seq.hpp"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <vector>

#include "romanov_m_matrix_ccs/common/include/common.hpp"

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
  const auto &a = GetInput().first;
  const auto &b = GetInput().second;
  auto &c = GetOutput();

  c.rows_num = a.rows_num;
  c.cols_num = b.cols_num;
  c.col_ptrs.assign(c.cols_num + 1, 0);

  std::vector<double> accumulator(a.rows_num, 0.0);
  std::vector<size_t> active_rows;
  std::vector<bool> row_mask(a.rows_num, false);

  for (size_t j = 0; j < b.cols_num; ++j) {
    c.col_ptrs[j] = c.vals.size();

    for (size_t kb = b.col_ptrs[j]; kb < b.col_ptrs[j + 1]; ++kb) {
      size_t k = b.row_inds[kb];
      double v_b = b.vals[kb];

      for (size_t ka = a.col_ptrs[k]; ka < a.col_ptrs[k + 1]; ++ka) {
        size_t i = a.row_inds[ka];
        if (!row_mask[i]) {
          row_mask[i] = true;
          active_rows.push_back(i);
        }
        accumulator[i] += a.vals[ka] * v_b;
      }
    }

    std::ranges::sort(active_rows);

    for (size_t row_idx : active_rows) {
      if (std::abs(accumulator[row_idx]) > 1e-12) {
        c.vals.push_back(accumulator[row_idx]);
        c.row_inds.push_back(row_idx);
      }
      accumulator[row_idx] = 0.0;
      row_mask[row_idx] = false;
    }
    active_rows.clear();
  }

  c.nnz = c.vals.size();
  c.col_ptrs[c.cols_num] = c.nnz;
  return true;
}

bool RomanovMMatrixCCSSeq::PostProcessingImpl() {
  return true;
}

}  // namespace romanov_m_matrix_ccs
