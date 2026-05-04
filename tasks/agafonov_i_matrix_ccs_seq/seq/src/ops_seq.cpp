#include "agafonov_i_matrix_ccs_seq/seq/include/ops_seq.hpp"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <utility>
#include <vector>

#include "agafonov_i_matrix_ccs_seq/common/include/common.hpp"

namespace agafonov_i_matrix_ccs_seq {

AgafonovIMatrixCCSSeq::AgafonovIMatrixCCSSeq(const InType &in) {
  SetTypeOfTask(GetStaticTypeOfTask());
  GetInput() = in;
}

bool AgafonovIMatrixCCSSeq::ValidationImpl() {
  const auto &left = GetInput().first;
  const auto &right = GetInput().second;

  return (left.cols_num == right.rows_num) && (left.col_ptrs.size() == left.cols_num + 1) &&
         (right.col_ptrs.size() == right.cols_num + 1);
}

bool AgafonovIMatrixCCSSeq::PreProcessingImpl() {
  GetOutput().vals.clear();
  GetOutput().row_inds.clear();
  return true;
}

bool AgafonovIMatrixCCSSeq::RunImpl() {
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
    c.col_ptrs[j] = static_cast<int>(c.vals.size());
    const auto b_col_start = static_cast<size_t>(b.col_ptrs[j]);
    const auto b_col_end = static_cast<size_t>(b.col_ptrs[j + 1]);

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
        c.vals.push_back(accumulator[row_idx]);
        c.row_inds.push_back(static_cast<int>(row_idx));
      }
      accumulator[row_idx] = 0.0;
      row_mask[row_idx] = false;
    }
    active_rows.clear();
  }

  c.nnz = c.vals.size();
  c.col_ptrs[c.cols_num] = static_cast<int>(c.nnz);
  return true;
}

bool AgafonovIMatrixCCSSeq::PostProcessingImpl() {
  return true;
}

}  // namespace agafonov_i_matrix_ccs_seq
