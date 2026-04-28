#include "guseva_crs/all/include/ops_all.hpp"

#include "guseva_crs/all/include/multiplier_all.hpp"
#include "guseva_crs/common/include/common.hpp"

namespace guseva_crs {

GusevaCRSMatMulAll::GusevaCRSMatMulAll(const InType &in) {
  SetTypeOfTask(GetStaticTypeOfTask());
  GetInput() = in;
  GetOutput();
}

bool GusevaCRSMatMulAll::ValidationImpl() {
  const auto &[a, b] = GetInput();
  return a.ncols == b.nrows;
}

bool GusevaCRSMatMulAll::PreProcessingImpl() {
  return true;
}

bool GusevaCRSMatMulAll::RunImpl() {
  const auto &[a, b] = GetInput();
  auto mult = MultiplierAll();
  GetOutput() = mult.Multiply(a, b);
  return true;
}

bool GusevaCRSMatMulAll::PostProcessingImpl() {
  return true;
}

}  // namespace guseva_crs
