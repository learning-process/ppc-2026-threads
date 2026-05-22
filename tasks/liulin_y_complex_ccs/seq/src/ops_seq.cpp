#include "liulin_y_complex_ccs/seq/include/ops_seq.hpp"

#include <algorithm>
#include <cmath>
#include <complex>
#include <vector>

#include "liulin_y_complex_ccs/common/include/common.hpp"

namespace liulin_y_complex_ccs {

LiulinYComplexCcs::LiulinYComplexCcs(const InType &in) : BaseTask() {
  SetTypeOfTask(GetStaticTypeOfTask());
  GetInput() = in;
}

bool LiulinYComplexCcs::ValidationImpl() {
  const auto &first = GetInput().first;
  const auto &second = GetInput().second;
  return first.count_cols == second.count_rows;
}

bool LiulinYComplexCcs::PreProcessingImpl() {
  const auto &first = GetInput().first;
  const auto &second = GetInput().second;

  auto &result = GetOutput();
  result.count_rows = first.count_rows;
  result.count_cols = second.count_cols;
  result.values.clear();
  result.row_index.clear();
  result.col_index.assign(result.count_cols + 1, 0);

  return true;
}

bool LiulinYComplexCcs::RunImpl() {
  const auto &first = GetInput().first;
  const auto &second = GetInput().second;
  auto &result = GetOutput();

  std::vector<std::complex<double>> dense_col(first.count_rows, {0.0, 0.0});
  std::vector<int> active_rows;
  std::vector<bool> is_active(first.count_rows, false);

  for (int j = 0; j < second.count_cols; ++j) {
    for (int kb = second.col_index[j]; kb < second.col_index[j + 1]; ++kb) {
      int k = second.row_index[kb];
      std::complex<double> b_val = second.values[kb];

      for (int ka = first.col_index[k]; ka < first.col_index[k + 1]; ++ka) {
        int i = first.row_index[ka];
        if (!is_active[i]) {
          is_active[i] = true;
          active_rows.push_back(i);
        }
        dense_col[i] += first.values[ka] * b_val;
      }
    }

    std::ranges::sort(active_rows);

    for (int i : active_rows) {
      if (std::abs(dense_col[i]) > 1e-15) {
        result.values.push_back(dense_col[i]);
        result.row_index.push_back(i);
      }
      dense_col[i] = {0.0, 0.0};
      is_active[i] = false;
    }
    result.col_index[j + 1] = static_cast<int>(result.values.size());
    active_rows.clear();
  }

  return true;
}

bool LiulinYComplexCcs::PostProcessingImpl() {
  return true;
}

}  // namespace liulin_y_complex_ccs
