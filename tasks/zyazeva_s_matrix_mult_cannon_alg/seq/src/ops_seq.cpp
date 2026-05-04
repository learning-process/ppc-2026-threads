#include "zyazeva_s_matrix_mult_cannon_alg/seq/include/ops_seq.hpp"

#include <cstddef>
#include <vector>

#include "zyazeva_s_matrix_mult_cannon_alg/common/include/common.hpp"

namespace zyazeva_s_matrix_mult_cannon_alg {

ZyazevaSMatrixMultCannonAlgSEQ::ZyazevaSMatrixMultCannonAlgSEQ(const InType &in) {
  SetTypeOfTask(GetStaticTypeOfTask());
  GetInput() = in;
  GetOutput() = {};
}

bool ZyazevaSMatrixMultCannonAlgSEQ::ValidationImpl() {
  const size_t sz = std::get<0>(GetInput());
  const auto &m1 = std::get<1>(GetInput());
  const auto &m2 = std::get<2>(GetInput());

  return sz > 0 && m1.size() == sz * sz && m2.size() == sz * sz;
}

bool ZyazevaSMatrixMultCannonAlgSEQ::PreProcessingImpl() {
  GetOutput() = {};
  return true;
}

bool ZyazevaSMatrixMultCannonAlgSEQ::RunImpl() {
  const auto sz = std::get<0>(GetInput());
  const auto &m1 = std::get<1>(GetInput());
  const auto &m2 = std::get<2>(GetInput());

  std::vector<double> res_m(sz * sz, 0.0);

  for (size_t row_idx = 0; row_idx < sz; ++row_idx) {
    for (size_t col_idx = 0; col_idx < sz; ++col_idx) {
      double accumulated = 0.0;
      for (size_t inner_idx = 0; inner_idx < sz; ++inner_idx) {
        accumulated += m1[(row_idx * sz) + inner_idx] * m2[(inner_idx * sz) + col_idx];
      }
      res_m[(row_idx * sz) + col_idx] = accumulated;
    }
  }

  GetOutput() = res_m;
  return true;
}

bool ZyazevaSMatrixMultCannonAlgSEQ::PostProcessingImpl() {
  return true;
}

}  // namespace zyazeva_s_matrix_mult_cannon_alg
