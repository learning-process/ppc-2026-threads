#include "golovanov_d_radix_merge/omp/include/ops_omp.hpp"

#include <vector>

#include "../include/radix_sort_omp.hpp"
#include "golovanov_d_radix_merge/common/include/common.hpp"

namespace golovanov_d_radix_merge {

GolovanovDRadixMergeOMP::GolovanovDRadixMergeOMP(const InType &in) {
  SetTypeOfTask(GetStaticTypeOfTask());
  GetInput() = in;
  GetOutput() = {};
}

bool GolovanovDRadixMergeOMP::ValidationImpl() {
  return !GetInput().empty();
}

bool GolovanovDRadixMergeOMP::PreProcessingImpl() {
  return true;
}

bool GolovanovDRadixMergeOMP::RunImpl() {
  std::vector<double> input = GetInput();
  RadixSortOMP::Sort(input);
  GetOutput() = input;
  return true;
}

bool GolovanovDRadixMergeOMP::PostProcessingImpl() {
  return true;
}

}  // namespace golovanov_d_radix_merge
