#include "mityaeva_radix/tbb/include/ops_tbb.hpp"

#include "mityaeva_radix/common/include/common.hpp"
#include "mityaeva_radix/tbb/include/sorter_tbb.hpp"

namespace mityaeva_radix {

MityaevaRadixTbb::MityaevaRadixTbb(const InType &in) {
  SetTypeOfTask(GetStaticTypeOfTask());
  GetInput() = in;
  GetOutput() = {};
}

bool MityaevaRadixTbb::ValidationImpl() {
  return !GetInput().empty();
}

bool MityaevaRadixTbb::PreProcessingImpl() {
  return true;
}

bool MityaevaRadixTbb::RunImpl() {
  auto &array = GetInput();
  SorterTbb::Sort(array);
  GetOutput() = array;
  return true;
}

bool MityaevaRadixTbb::PostProcessingImpl() {
  return true;
}

}  // namespace mityaeva_radix
