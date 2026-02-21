#include "kulik_a_mat_mul_double_ccs/seq/include/ops_seq.hpp"

#include <algorithm>
#include <cstddef>
#include <tuple>
#include <vector>

#include "kulik_a_mat_mul_double_ccs/common/include/common.hpp"

namespace kulik_a_mat_mul_double_ccs {

KulikAMatMulDoubleCcsSEQ::KulikAMatMulDoubleCcsSEQ(const InType &in) {
  SetTypeOfTask(GetStaticTypeOfTask());
  GetInput() = in;
}

bool KulikAMatMulDoubleCcsSEQ::ValidationImpl() {
  const auto &a = std::get<0>(GetInput());
  const auto &b = std::get<1>(GetInput());
  return (a.m == b.n);
}

bool KulikAMatMulDoubleCcsSEQ::PreProcessingImpl() {
  return true;
}

bool KulikAMatMulDoubleCcsSEQ::RunImpl() {
  const auto &a = std::get<0>(GetInput());
  const auto &b = std::get<1>(GetInput());
  OutType &c = GetOutput();
  c.value.clear();
  c.row.clear();
  c.n = a.n;
  c.m = b.m;
  c.col_ind.assign(c.m + 1, 0);                //
  std::vector<double> accum(a.n, 0.0);         //
  std::vector<bool> nz_elem_rows(a.n, false);  //
  std::vector<size_t> nnz_rows;
  for (size_t j = 0; j < b.m; ++j) {
    c.col_ind[j] = c.value.size();
    for (size_t k = b.col_ind[j]; k < b.col_ind[j + 1]; ++k) {
      size_t ind = b.row[k];
      double b_val = b.value[k];
      for (size_t zc = a.col_ind[ind]; zc < a.col_ind[ind + 1]; ++zc) {
        size_t i = a.row[zc];
        double a_val = a.value[zc];
        accum[i] += a_val * b_val;
        if (!nz_elem_rows[i]) {
          nz_elem_rows[i] = true;
          nnz_rows.push_back(i);
        }
      }
    }
    std::ranges::sort(nnz_rows);
    ;
    for (size_t i : nnz_rows) {
      if (accum[i] != 0.0) {
        c.row.push_back(i);
        c.value.push_back(accum[i]);
      }
      accum[i] = 0.0;
      nz_elem_rows[i] = false;
    }
    nnz_rows.clear();
  }
  c.nz = c.value.size();
  c.col_ind[b.m] = c.nz;

  return true;
}

bool KulikAMatMulDoubleCcsSEQ::PostProcessingImpl() {
  return true;
}

}  // namespace kulik_a_mat_mul_double_ccs
