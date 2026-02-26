#include "tsyplakov_k_mul_double_crs_matrix/tbb/include/ops_tbb.hpp"

#include <tbb/tbb.h>

#include <atomic>
#include <numeric>
#include <util/include/util.hpp>
#include <vector>

#include "oneapi/tbb/parallel_for.h"
#include "tsyplakov_k_mul_double_crs_matrix/common/include/common.hpp"

namespace tsyplakov_k_mul_double_crs_matrix {

TsyplakovMulCrsMatrixTBB::TsyplakovMulCrsMatrixTBB(const InType &in) {
  SetTypeOfTask(GetStaticTypeOfTask());
  GetInput() = in;
}

bool TsyplakovMulCrsMatrixTBB::ValidationImpl() {
  return true;
}

bool TsyplakovMulCrsMatrixTBB::PreProcessingImpl() {
  return true;
}

bool TsyplakovMulCrsMatrixTBB::RunImpl() {
  return true;
}

bool TsyplakovMulCrsMatrixTBB::PostProcessingImpl() {
  return true;
}

}  // namespace tsyplakov_k_mul_double_crs_matrix
