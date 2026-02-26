#include "liulin_y_complex_ccs/seq/include/ops_seq.hpp"

#include <algorithm>
#include <complex>
#include <numeric>
#include <vector>


#include "liulin_y_complex_ccs/common/include/common.hpp"
#include "util/include/util.hpp"

namespace liulin_y_complex_ccs {

LiulinYComplexCcs::LiulinYComplexCcs(const InType &in) : BaseTask() {
  SetTypeOfTask(GetStaticTypeOfTask());
  GetInput() = in;
}

bool LiulinYComplexCcs::ValidationImpl() {
  const auto& A = GetInput().first;
  const auto& B = GetInput().second;
  return A.count_cols == B.count_rows;
}

bool LiulinYComplexCcs::PreProcessingImpl() {
  const auto& A = GetInput().first;
  const auto& B = GetInput().second;
  
  auto& Result = GetOutput();
  Result.count_rows = A.count_rows;
  Result.count_cols = B.count_cols;
  Result.values.clear();
  Result.row_index.clear();
  Result.col_index.assign(Result.count_cols + 1, 0);
  
  return true;
}

bool LiulinYComplexCcs::RunImpl() {
  const auto& A = GetInput().first;
  const auto& B = GetInput().second;
  auto& Result = GetOutput();

  std::vector<std::complex<double>> dense_col(A.count_rows, {0.0, 0.0});
  std::vector<int> active_rows;
  std::vector<bool> is_active(A.count_rows, false);

  for (int j = 0; j < B.count_cols; ++j) {
    for (int kb = B.col_index[j]; kb < B.col_index[j + 1]; ++kb) {
      int k = B.row_index[kb];              
      std::complex<double> b_val = B.values[kb];

      for (int ka = A.col_index[k]; ka < A.col_index[k + 1]; ++ka) {
        int i = A.row_index[ka];        
        if (!is_active[i]) {
          is_active[i] = true;
          active_rows.push_back(i);
        }
        dense_col[i] += A.values[ka] * b_val;
      }
    }

    std::sort(active_rows.begin(), active_rows.end());

    for (int i : active_rows) {
      if (std::abs(dense_col[i]) > 1e-15) { 
        Result.values.push_back(dense_col[i]);
        Result.row_index.push_back(i);
      }
      dense_col[i] = {0.0, 0.0};
      is_active[i] = false;
    }
    Result.col_index[j + 1] = static_cast<int>(Result.values.size());
    active_rows.clear();
  }

  return true;
}

bool LiulinYComplexCcs::PostProcessingImpl() {
  return true;
}

}  // namespace liulin_y_complex_ccs