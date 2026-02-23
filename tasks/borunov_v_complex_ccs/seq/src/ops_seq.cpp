#include "borunov_v_complex_ccs/seq/include/ops_seq.hpp"

#include <algorithm>
#include <complex>
#include <vector>

#include "borunov_v_complex_ccs/common/include/common.hpp"
#include "util/include/util.hpp"

namespace borunov_v_complex_ccs_seq {

BorunovVComplexCcsSEQ::BorunovVComplexCcsSEQ(const InType &in) {
  SetTypeOfTask(GetStaticTypeOfTask());
  GetInput() = in;
  GetOutput().resize(1);
}

bool BorunovVComplexCcsSEQ::ValidationImpl() {
  if (GetInput().size() != 2) {
    return false;
  }
  const auto& A = GetInput()[0];
  const auto& B = GetInput()[1];
  if (A.num_cols != B.num_rows) {
    return false;
  }
  if (A.col_ptrs.size() != static_cast<size_t>(A.num_cols + 1) ||
      B.col_ptrs.size() != static_cast<size_t>(B.num_cols + 1)) {
    return false;
  }
  return true;
}

bool BorunovVComplexCcsSEQ::PreProcessingImpl() {
  const auto& A = GetInput()[0];
  const auto& B = GetInput()[1];
  auto& C = GetOutput()[0];

  C.num_rows = A.num_rows;
  C.num_cols = B.num_cols;
  C.col_ptrs.assign(C.num_cols + 1, 0);
  C.values.clear();
  C.row_indices.clear();

  return true;
}

bool BorunovVComplexCcsSEQ::RunImpl() {
  const auto& A = GetInput()[0];
  const auto& B = GetInput()[1];
  auto& C = GetOutput()[0];

  std::vector<std::complex<double>> col_accumulator(A.num_rows, {0.0, 0.0});
  std::vector<int> non_zero_indices;
  std::vector<bool> is_non_zero(A.num_rows, false);

  for (int j = 0; j < B.num_cols; ++j) {
    for (int b_idx = B.col_ptrs[j]; b_idx < B.col_ptrs[j + 1]; ++b_idx) {
      int p = B.row_indices[b_idx];
      std::complex<double> b_val = B.values[b_idx];

      for (int a_idx = A.col_ptrs[p]; a_idx < A.col_ptrs[p + 1]; ++a_idx) {
        int i = A.row_indices[a_idx];
        std::complex<double> a_val = A.values[a_idx];

        if (!is_non_zero[i]) {
          is_non_zero[i] = true;
          non_zero_indices.push_back(i);
        }
        col_accumulator[i] += a_val * b_val;
      }
    }

    std::sort(non_zero_indices.begin(), non_zero_indices.end());

    for (int i : non_zero_indices) {
      if (std::abs(col_accumulator[i]) > 1e-9) {
        C.values.push_back(col_accumulator[i]);
        C.row_indices.push_back(i);
      }
      col_accumulator[i] = {0.0, 0.0};
      is_non_zero[i] = false;
    }
    non_zero_indices.clear();

    C.col_ptrs[j + 1] = C.values.size();
  }

  return true;
}

bool BorunovVComplexCcsSEQ::PostProcessingImpl() {
  return true;
}

}  // namespace borunov_v_complex_ccs_seq
