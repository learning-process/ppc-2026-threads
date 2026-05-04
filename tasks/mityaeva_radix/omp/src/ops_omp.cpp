#include "mityaeva_radix/omp/include/ops_omp.hpp"

#include "mityaeva_radix/common/include/common.hpp"
#include "mityaeva_radix/omp/include/sorter_omp.hpp"

namespace mityaeva_radix {

MityaevaRadixOmp::MityaevaRadixOmp(const InType &in) {
  SetTypeOfTask(GetStaticTypeOfTask());
  GetInput() = in;
  GetOutput() = {};
}

bool MityaevaRadixOmp::ValidationImpl() {
  return !GetInput().empty();
}

bool MityaevaRadixOmp::PreProcessingImpl() {
  return true;
}

bool MityaevaRadixOmp::RunImpl() {
  auto &array = GetInput();
  SorterOmp::Sort(array);
  GetOutput() = array;
  return true;
}

bool MityaevaRadixOmp::PostProcessingImpl() {
  return true;
}

}  // namespace mityaeva_radix
