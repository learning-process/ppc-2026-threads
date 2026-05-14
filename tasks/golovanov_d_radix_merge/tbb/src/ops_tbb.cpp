#include "golovanov_d_radix_merge/tbb/include/ops_tbb.hpp"

#include <vector>

#include "../include/radix_sort_tbb.hpp"
#include "golovanov_d_radix_merge/common/include/common.hpp"

namespace golovanov_d_radix_merge {

GolovanovDRadixMergeTBB::GolovanovDRadixMergeTBB(const InType &in) {
  SetTypeOfTask(GetStaticTypeOfTask());
  GetInput() = in;
  GetOutput() = {};
}

bool GolovanovDRadixMergeTBB::ValidationImpl() {
  return !GetInput().empty();
}

bool GolovanovDRadixMergeTBB::PreProcessingImpl() {
  return true;
}

bool GolovanovDRadixMergeTBB::RunImpl() {
  std::vector<double> input = GetInput();
  RadixSortTBB::Sort(input);
  GetOutput() = input;
  return true;
}

bool GolovanovDRadixMergeTBB::PostProcessingImpl() {
  return true;
}

}  // namespace golovanov_d_radix_merge
