#include "guseva_crs/tbb/include/ops_tbb.hpp"

#include "guseva_crs/common/include/common.hpp"
#include "guseva_crs/tbb/include/multiplier_tbb.hpp"

namespace guseva_crs {

GusevaCRSMatMulTbb::GusevaCRSMatMulTbb(const InType &in) {
  SetTypeOfTask(GetStaticTypeOfTask());
  GetInput() = in;
  GetOutput();
}

bool GusevaCRSMatMulTbb::ValidationImpl() {
  const auto &[a, b] = GetInput();
  return a.ncols == b.nrows;
}

bool GusevaCRSMatMulTbb::PreProcessingImpl() {
  return true;
}

bool GusevaCRSMatMulTbb::RunImpl() {
  const auto &[a, b] = GetInput();
  auto mult = MultiplierTbb();
  GetOutput() = mult.Multiply(a, b);
  return true;
}

bool GusevaCRSMatMulTbb::PostProcessingImpl() {
  return true;
}

}  // namespace guseva_crs
