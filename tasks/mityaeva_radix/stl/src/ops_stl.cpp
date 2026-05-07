#include "mityaeva_radix/stl/include/ops_stl.hpp"

#include "mityaeva_radix/common/include/common.hpp"
#include "mityaeva_radix/stl/include/sorter_stl.hpp"

namespace mityaeva_radix {

MityaevaRadixStl::MityaevaRadixStl(const InType &in) {
  SetTypeOfTask(GetStaticTypeOfTask());
  GetInput() = in;
  GetOutput() = {};
}

bool MityaevaRadixStl::ValidationImpl() {
  return !GetInput().empty();
}

bool MityaevaRadixStl::PreProcessingImpl() {
  return true;
}

bool MityaevaRadixStl::RunImpl() {
  auto &array = GetInput();
  SorterStl::Sort(array);
  GetOutput() = array;
  return true;
}

bool MityaevaRadixStl::PostProcessingImpl() {
  return true;
}

}  // namespace mityaeva_radix
