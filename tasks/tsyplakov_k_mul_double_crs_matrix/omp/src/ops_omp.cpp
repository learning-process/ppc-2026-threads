#include "tsyplakov_k_mul_double_crs_matrix/omp/include/ops_omp.hpp"

#include <atomic>
#include <numeric>
#include <vector>

#include "tsyplakov_k_mul_double_crs_matrix/common/include/common.hpp"
#include "util/include/util.hpp"

namespace tsyplakov_k_mul_double_crs_matrix {

TsyplakovMulCrsMatrixOMP::TsyplakovMulCrsMatrixOMP(const InType &in) {
  SetTypeOfTask(GetStaticTypeOfTask());
  GetInput() = in;
}

bool TsyplakovMulCrsMatrixOMP::ValidationImpl() {
  return true;
}

bool TsyplakovMulCrsMatrixOMP::PreProcessingImpl() {
  return true;
}

bool TsyplakovMulCrsMatrixOMP::RunImpl() {
  return true;
}

bool TsyplakovMulCrsMatrixOMP::PostProcessingImpl() {
  return true;
}

}  // namespace tsyplakov_k_mul_double_crs_matrix
