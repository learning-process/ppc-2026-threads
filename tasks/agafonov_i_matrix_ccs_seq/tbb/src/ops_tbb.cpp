#include "agafonov_i_matrix_ccs_seq/tbb/include/ops_tbb.hpp"

#include <tbb/blocked_range.h>
#include <tbb/parallel_for.h>

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <vector>

#include "agafonov_i_matrix_ccs_seq/common/include/common.hpp"

namespace agafonov_i_matrix_ccs_seq {

AgafonovIMatrixCCSTBB::AgafonovIMatrixCCSTBB(const InType &in) {
  SetTypeOfTask(GetStaticTypeOfTask());
  GetInput() = in;
}

bool AgafonovIMatrixCCSTBB::ValidationImpl() {
  const auto &left = GetInput().first;
  const auto &right = GetInput().second;
  return (left.cols_num == right.rows_num) && (left.col_ptrs.size() == left.cols_num + 1) &&
         (right.col_ptrs.size() == right.cols_num + 1);
}

bool AgafonovIMatrixCCSTBB::PreProcessingImpl() {
  GetOutput().vals.clear();
  GetOutput().row_inds.clear();
  return true;
}

void AgafonovIMatrixCCSTBB::ProcessColumn(size_t j, const InType::first_type &a, const InType::second_type &b,
                                          std::vector<double> &accumulator, std::vector<size_t> &active_rows,
                                          std::vector<bool> &row_mask, std::vector<double> &local_v,
                                          std::vector<int> &local_r) {
  const auto b_col_start = static_cast<size_t>(b.col_ptrs[j]);
  const auto b_col_end = static_cast<size_t>(b.col_ptrs[j + 1]);
  if (b_col_start == b_col_end) {
    return;
  }

  for (size_t kb = b_col_start; kb < b_col_end; ++kb) {
    const auto k = static_cast<size_t>(b.row_inds[kb]);
    const double v_b = b.vals[kb];
    const auto a_col_start = static_cast<size_t>(a.col_ptrs[k]);
    const auto a_col_end = static_cast<size_t>(a.col_ptrs[k + 1]);

    for (size_t ka = a_col_start; ka < a_col_end; ++ka) {
      const auto i = static_cast<size_t>(a.row_inds[ka]);
      if (!row_mask[i]) {
        row_mask[i] = true;
        active_rows.push_back(i);
      }
      accumulator[i] += a.vals[ka] * v_b;
    }
  }

  std::ranges::sort(active_rows);

  for (const auto row_idx : active_rows) {
    if (std::abs(accumulator[row_idx]) > 1e-15) {
      local_v.push_back(accumulator[row_idx]);
      local_r.push_back(static_cast<int>(row_idx));
    }
    row_mask[row_idx] = false;
    accumulator[row_idx] = 0.0;
  }
  active_rows.clear();
}

bool AgafonovIMatrixCCSTBB::RunImpl() {
  const auto &a = GetInput().first;
  const auto &b = GetInput().second;
  auto &c = GetOutput();

  c.rows_num = a.rows_num;
  c.cols_num = b.cols_num;
  c.col_ptrs.assign(c.cols_num + 1, 0);

  std::vector<std::vector<double>> local_vals(b.cols_num);
  std::vector<std::vector<int>> local_rows(b.cols_num);

  tbb::parallel_for(tbb::blocked_range<size_t>(0, b.cols_num), [&](const tbb::blocked_range<size_t> &range) {
    std::vector<double> accumulator(a.rows_num, 0.0);
    std::vector<size_t> active_rows;
    std::vector<bool> row_mask(a.rows_num, false);

    for (size_t j = range.begin(); j < range.end(); ++j) {
      ProcessColumn(j, a, b, accumulator, active_rows, row_mask, local_vals[j], local_rows[j]);
    }
  });

  int current_nnz = 0;
  for (size_t j = 0; j < b.cols_num; ++j) {
    c.col_ptrs[j] = current_nnz;
    c.vals.insert(c.vals.end(), local_vals[j].begin(), local_vals[j].end());
    c.row_inds.insert(c.row_inds.end(), local_rows[j].begin(), local_rows[j].end());
    current_nnz += static_cast<int>(local_vals[j].size());
  }
  c.col_ptrs[b.cols_num] = current_nnz;
  c.nnz = current_nnz;
  return true;
}

bool AgafonovIMatrixCCSTBB::PostProcessingImpl() {
  return true;
}

}  // namespace agafonov_i_matrix_ccs_seq
