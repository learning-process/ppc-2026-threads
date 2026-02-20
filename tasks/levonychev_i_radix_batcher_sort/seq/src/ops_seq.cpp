#include "levonychev_i_radix_batcher_sort/seq/include/ops_seq.hpp"

#include <numeric>
#include <vector>

#include "levonychev_i_radix_batcher_sort/common/include/common.hpp"
#include "util/include/util.hpp"

namespace levonychev_i_radix_batcher_sort {

LevonychevIRadixBatcherSortSEQ::LevonychevIRadixBatcherSortSEQ(const InType &in) {
  SetTypeOfTask(GetStaticTypeOfTask());
  GetInput() = in;
  GetOutput() = 0;
}

bool LevonychevIRadixBatcherSortSEQ::ValidationImpl() {
  return GetInput() >= 0 || GetInput() < 0;
}

bool LevonychevIRadixBatcherSortSEQ::PreProcessingImpl() {
  return GetOutput() == 0;
}

bool LevonychevIRadixBatcherSortSEQ::RunImpl() {
  return GetOutput() == 0;
}

bool LevonychevIRadixBatcherSortSEQ::PostProcessingImpl() {
  return GetOutput() == 0;
}

}  // namespace levonychev_i_radix_batcher_sort
