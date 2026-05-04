#include "guseva_crs/seq/include/ops_seq.hpp"

#include "guseva_crs/common/include/common.hpp"
#include "guseva_crs/seq/include/multiplier_seq.hpp"

namespace guseva_crs {

GusevaCRSMatMulSeq::GusevaCRSMatMulSeq(const InType &in) {
  SetTypeOfTask(GetStaticTypeOfTask());
  GetInput() = in;
  GetOutput();
}

bool GusevaCRSMatMulSeq::ValidationImpl() {
  const auto &[a, b] = GetInput();
  return a.ncols == b.nrows;
}

bool GusevaCRSMatMulSeq::PreProcessingImpl() {
  return true;
}

bool GusevaCRSMatMulSeq::RunImpl() {
  const auto &[a, b] = GetInput();
  auto mult = MultiplierSeq();
  GetOutput() = mult.Multiply(a, b);
  return true;
}

bool GusevaCRSMatMulSeq::PostProcessingImpl() {
  return true;
}

}  // namespace guseva_crs
