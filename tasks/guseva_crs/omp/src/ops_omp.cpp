#include "guseva_crs/omp/include/ops_omp.hpp"

#include "guseva_crs/common/include/common.hpp"
#include "guseva_crs/omp/include/multiplier_omp.hpp"

namespace guseva_crs {

GusevaCRSMatMulOmp::GusevaCRSMatMulOmp(const InType &in) {
  SetTypeOfTask(GetStaticTypeOfTask());
  GetInput() = in;
  GetOutput();
}

bool GusevaCRSMatMulOmp::ValidationImpl() {
  const auto &[a, b] = GetInput();
  return a.ncols == b.nrows;
}

bool GusevaCRSMatMulOmp::PreProcessingImpl() {
  return true;
}

bool GusevaCRSMatMulOmp::RunImpl() {
  const auto &[a, b] = GetInput();
  auto mult = MultiplierOmp();
  GetOutput() = mult.Multiply(a, b);
  return true;
}

bool GusevaCRSMatMulOmp::PostProcessingImpl() {
  return true;
}

}  // namespace guseva_crs
