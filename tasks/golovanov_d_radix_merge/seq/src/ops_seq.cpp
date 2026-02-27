#include "golovanov_d_radix_merge/seq/include/ops_seq.hpp"

#include <numeric>
#include <vector>

#include "golovanov_d_radix_merge/common/include/common.hpp"
#include "util/include/util.hpp"

namespace golovanov_d_radix_merge {

GolovanovDRadixMergeSEQ::GolovanovDRadixMergeSEQ(const InType &in) {
  SetTypeOfTask(GetStaticTypeOfTask());
  GetInput() = in;
  GetOutput() = 0;
}

bool GolovanovDRadixMergeSEQ::ValidationImpl() {
  return (GetInput() > 0) && (GetOutput() == 0);
}

bool GolovanovDRadixMergeSEQ::PreProcessingImpl() {
  return true;
}

bool GolovanovDRadixMergeSEQ::RunImpl() {
  
}

bool GolovanovDRadixMergeSEQ::PostProcessingImpl() {
  return true;
}

}  // namespace golovanov_d_radix_merge
