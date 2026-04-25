#include "golovanov_d_radix_merge/stl/include/ops_stl.hpp"

#include <vector>

#include "../include/radix_sort_stl.hpp"
#include "golovanov_d_radix_merge/common/include/common.hpp"

namespace golovanov_d_radix_merge {

GolovanovDRadixMergeSTL::GolovanovDRadixMergeSTL(const InType &in) {
  SetTypeOfTask(GetStaticTypeOfTask());
  GetInput() = in;
  GetOutput() = {};
}

bool GolovanovDRadixMergeSTL::ValidationImpl() {
  return !GetInput().empty();
}

bool GolovanovDRadixMergeSTL::PreProcessingImpl() {
  return true;
}

bool GolovanovDRadixMergeSTL::RunImpl() {
  std::vector<double> input = GetInput();
  RadixSortSTL::Sort(input);
  GetOutput() = input;
  return true;
}

bool GolovanovDRadixMergeSTL::PostProcessingImpl() {
  return true;
}

}  // namespace golovanov_d_radix_merge
