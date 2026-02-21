#include "example_threads/seq/include/ops_seq.hpp"

#include <numeric>
#include <vector>
#include <tuple>

#include "example_threads/common/include/common.hpp"
#include "util/include/util.hpp"

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
  c.n = a.n;
  c.m = b.m;
  c.col_ind.resize(c.m + 1, 0);


  return true;
}

bool KulikAMatMulDoubleCcsSEQ::PostProcessingImpl() {
  GetOutput() -= GetInput();
  return GetOutput() > 0;
}

}  // namespace kulik_a_mat_mul_double_ccs
