#include "golovanov_d_radix_merge/seq/include/ops_seq.hpp"

#include <numeric>
#include <vector>

#include "../include/radix_sort.hpp"
#include "golovanov_d_radix_merge/common/include/common.hpp"
#include "util/include/util.hpp"

namespace golovanov_d_radix_merge {

GolovanovDRadixMergeSEQ::GolovanovDRadixMergeSEQ(const InType &in) {
  SetTypeOfTask(GetStaticTypeOfTask());
  GetInput() = in;
  GetOutput() = {};
}

bool GolovanovDRadixMergeSEQ::ValidationImpl() {
  return !GetInput().empty();
}

bool GolovanovDRadixMergeSEQ::PreProcessingImpl() {
  return true;
}

bool GolovanovDRadixMergeSEQ::RunImpl() {
  std::vector<double> input = GetInput();
  RadixSort::Sort(input);
  GetOutput() = input;
  return true;
}

bool GolovanovDRadixMergeSEQ::PostProcessingImpl() {
  return true;
}

}  // namespace golovanov_d_radix_merge
