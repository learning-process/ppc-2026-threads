#include "borunov_v_complex_ccs/seq/include/ops_seq.hpp"

#include <algorithm>
#include <cmath>
#include <complex>
#include <cstddef>
#include <vector>

#include "borunov_v_complex_ccs/common/include/common.hpp"

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
  const auto &a = GetInput()[0];
  const auto &b = GetInput()[1];
  if (a.num_cols != b.num_rows) {
    return false;
  }
  if (a.col_ptrs.size() != static_cast<size_t>(a.num_cols) + 1 ||
      b.col_ptrs.size() != static_cast<size_t>(b.num_cols) + 1) {
    return false;
  }
  return true;
}

bool BorunovVComplexCcsSEQ::PreProcessingImpl() {
  const auto &a = GetInput()[0];
  const auto &b = GetInput()[1];
  auto &c = GetOutput()[0];

  c.num_rows = a.num_rows;
  c.num_cols = b.num_cols;
  c.col_ptrs.assign(c.num_cols + 1, 0);
  c.values.clear();
  c.row_indices.clear();

  return true;
}

bool BorunovVComplexCcsSEQ::RunImpl() {
  const auto &a = GetInput()[0];
  const auto &b = GetInput()[1];
  auto &c = GetOutput()[0];

  std::vector<std::complex<double>> col_accumulator(a.num_rows, {0.0, 0.0});
  std::vector<int> non_zero_indices;
  std::vector<bool> is_non_zero(a.num_rows, false);

  for (int j = 0; j < b.num_cols; ++j) {
    for (int b_idx = b.col_ptrs[j]; b_idx < b.col_ptrs[j + 1]; ++b_idx) {
      int p = b.row_indices[b_idx];
      std::complex<double> b_val = b.values[b_idx];

      for (int a_idx = a.col_ptrs[p]; a_idx < a.col_ptrs[p + 1]; ++a_idx) {
        int i = a.row_indices[a_idx];
        std::complex<double> a_val = a.values[a_idx];

        if (!is_non_zero[i]) {
          is_non_zero[i] = true;
          non_zero_indices.push_back(i);
        }
        col_accumulator[i] += a_val * b_val;
      }
    }

    std::ranges::sort(non_zero_indices);

    for (int i : non_zero_indices) {
      if (std::abs(col_accumulator[i]) > 1e-9) {
        c.values.push_back(col_accumulator[i]);
        c.row_indices.push_back(i);
      }
      col_accumulator[i] = {0.0, 0.0};
      is_non_zero[i] = false;
    }
    non_zero_indices.clear();

    c.col_ptrs[j + 1] = static_cast<int>(c.values.size());
  }

  return true;
}

bool BorunovVComplexCcsSEQ::PostProcessingImpl() {
  return true;
}

}  // namespace borunov_v_complex_ccs_seq
