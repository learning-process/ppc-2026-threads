#include "tsyplakov_k_mul_double_crs_matrix/stl/include/ops_stl.hpp"

#include <atomic>
#include <numeric>
#include <thread>
#include <vector>

#include "tsyplakov_k_mul_double_crs_matrix/common/include/common.hpp"
#include "util/include/util.hpp"

namespace tsyplakov_k_mul_double_crs_matrix {

TsyplakovMulCrsMatrixSTL::TsyplakovMulCrsMatrixSTL(const InType &in) {
  SetTypeOfTask(GetStaticTypeOfTask());
  GetInput() = in;
}

bool TsyplakovMulCrsMatrixSTL::ValidationImpl() {
  return true;
}

bool TsyplakovMulCrsMatrixSTL::PreProcessingImpl() {
  return true;
}

bool TsyplakovMulCrsMatrixSTL::RunImpl() {
  return true;
}

bool TsyplakovMulCrsMatrixSTL::PostProcessingImpl() {
  return true;
}

}  // namespace tsyplakov_k_mul_double_crs_matrix
