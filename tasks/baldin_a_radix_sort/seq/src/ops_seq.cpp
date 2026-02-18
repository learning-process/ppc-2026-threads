#include "baldin_a_radix_sort/seq/include/ops_seq.hpp"

#include <numeric>
#include <vector>

#include "baldin_a_radix_sort/common/include/common.hpp"
#include "util/include/util.hpp"

namespace baldin_a_radix_sort {

BaldinARadixSortSEQ::BaldinARadixSortSEQ(const InType &in) {
  SetTypeOfTask(GetStaticTypeOfTask());
  GetInput() = in;
  GetOutput() = {};
}

bool BaldinARadixSortSEQ::ValidationImpl() {
  return true;
}

bool BaldinARadixSortSEQ::PreProcessingImpl() {
  return true;
}

bool BaldinARadixSortSEQ::RunImpl() {
  
}

bool BaldinARadixSortSEQ::PostProcessingImpl() {
  return true;
}

}  // namespace baldin_a_radix_sort
