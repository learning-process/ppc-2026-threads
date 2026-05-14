#include "mityaeva_radix/seq/include/ops_seq.hpp"

#include "mityaeva_radix/common/include/common.hpp"
#include "mityaeva_radix/seq/include/sorter_seq.hpp"

namespace mityaeva_radix {

MityaevaRadixSeq::MityaevaRadixSeq(const InType &in) {
  SetTypeOfTask(GetStaticTypeOfTask());
  GetInput() = in;
  GetOutput() = {};
}

bool MityaevaRadixSeq::ValidationImpl() {
  return !GetInput().empty();
}

bool MityaevaRadixSeq::PreProcessingImpl() {
  return true;
}

bool MityaevaRadixSeq::RunImpl() {
  auto &array = GetInput();
  SorterSeq::LSDSortDouble(array);
  GetOutput() = array;
  return true;
}

bool MityaevaRadixSeq::PostProcessingImpl() {
  return true;
}

}  // namespace mityaeva_radix
