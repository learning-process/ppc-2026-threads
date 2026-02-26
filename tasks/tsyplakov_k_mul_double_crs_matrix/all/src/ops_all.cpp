#include "tsyplakov_k_mul_double_crs_matrix/all/include/ops_all.hpp"

#include <mpi.h>

#include <atomic>
#include <numeric>
#include <thread>
#include <vector>

#include "oneapi/tbb/parallel_for.h"
#include "tsyplakov_k_mul_double_crs_matrix/common/include/common.hpp"
#include "util/include/util.hpp"

namespace tsyplakov_k_mul_double_crs_matrix {

TsyplakovKMulCrsMatrixALL::TsyplakovKMulCrsMatrixALL(const InType &in) {
  SetTypeOfTask(GetStaticTypeOfTask());
  GetInput() = in;
}

bool TsyplakovKMulCrsMatrixALL::ValidationImpl() {
  return true;
}

bool TsyplakovKMulCrsMatrixALL::PreProcessingImpl() {
  return true;
}

bool TsyplakovKMulCrsMatrixALL::RunImpl() {
  return true;
}

bool TsyplakovKMulCrsMatrixALL::PostProcessingImpl() {
  return true;
}

}  // namespace tsyplakov_k_mul_double_crs_matrix
