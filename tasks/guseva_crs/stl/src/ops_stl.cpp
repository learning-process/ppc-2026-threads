#include "guseva_crs/stl/include/ops_stl.hpp"

#include "guseva_crs/common/include/common.hpp"
#include "guseva_crs/stl/include/multiplier_stl.hpp"

namespace guseva_crs {

GusevaCRSMatMulStl::GusevaCRSMatMulStl(const InType &in) {
  SetTypeOfTask(GetStaticTypeOfTask());
  GetInput() = in;
  GetOutput();
}

bool GusevaCRSMatMulStl::ValidationImpl() {
  const auto &[a, b] = GetInput();
  return a.ncols == b.nrows;
}

bool GusevaCRSMatMulStl::PreProcessingImpl() {
  return true;
}

bool GusevaCRSMatMulStl::RunImpl() {
  const auto &[a, b] = GetInput();
  auto mult = MultiplierStl();
  GetOutput() = mult.Multiply(a, b);
  return true;
}

bool GusevaCRSMatMulStl::PostProcessingImpl() {
  return true;
}

}  // namespace guseva_crs
