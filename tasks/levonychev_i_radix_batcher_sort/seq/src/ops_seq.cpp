#include "levonychev_i_radix_batcher_sort/seq/include/ops_seq.hpp"

#include <numeric>
#include <vector>

#include "levonychev_i_radix_batcher_sort/common/include/common.hpp"
#include "util/include/util.hpp"

namespace levonychev_i_radix_batcher_sort {

LevonychevIRadixBatcherSortSEQ::LevonychevIRadixBatcherSortSEQ(const InType &in) {
  SetTypeOfTask(GetStaticTypeOfTask());
  GetInput() = in;
}

bool LevonychevIRadixBatcherSortSEQ::ValidationImpl() {
  return true;
}

bool LevonychevIRadixBatcherSortSEQ::PreProcessingImpl() {
  return true;
}

bool LevonychevIRadixBatcherSortSEQ::RunImpl() {
  
  GetOutput() = GetInput();
  std::sort(GetOutput().begin(), GetOutput().end());
  for (const auto& i : GetOutput())
    std::cout << i << ' ';
  std::cout << std::endl;
  return true;
}

bool LevonychevIRadixBatcherSortSEQ::PostProcessingImpl() {
  return true;
}

}  // namespace levonychev_i_radix_batcher_sort
