#include "dorofeev_i_bitwise_sort_double_eo_batcher_merge/tbb/include/ops_tbb.hpp"

#include <tbb/tbb.h>

#include <atomic>
#include <numeric>
#include <util/include/util.hpp>
#include <vector>

#include "dorofeev_i_bitwise_sort_double_eo_batcher_merge/common/include/common.hpp"
#include "oneapi/tbb/parallel_for.h"

namespace dorofeev_i_bitwise_sort_double_eo_batcher_merge {

DorofeevIBitwiseSortDoubleEOBatcherMergeTBB::DorofeevIBitwiseSortDoubleEOBatcherMergeTBB(const InType &in) {
  SetTypeOfTask(GetStaticTypeOfTask());
  GetInput() = in;
  GetOutput() = 0;
}

bool DorofeevIBitwiseSortDoubleEOBatcherMergeTBB::ValidationImpl() {
  return (GetInput() > 0) && (GetOutput() == 0);
}

bool DorofeevIBitwiseSortDoubleEOBatcherMergeTBB::PreProcessingImpl() {
  GetOutput() = 2 * GetInput();
  return GetOutput() > 0;
}

bool DorofeevIBitwiseSortDoubleEOBatcherMergeTBB::RunImpl() {
  for (InType i = 0; i < GetInput(); i++) {
    for (InType j = 0; j < GetInput(); j++) {
      for (InType k = 0; k < GetInput(); k++) {
        std::vector<InType> tmp(i + j + k, 1);
        GetOutput() += std::accumulate(tmp.begin(), tmp.end(), 0);
        GetOutput() -= i + j + k;
      }
    }
  }

  const int num_threads = ppc::util::GetNumThreads();
  GetOutput() *= num_threads;

  std::atomic<int> counter(0);
  tbb::parallel_for(0, ppc::util::GetNumThreads(), [&](int /*i*/) { counter++; });

  GetOutput() /= counter;
  return GetOutput() > 0;
}

bool DorofeevIBitwiseSortDoubleEOBatcherMergeTBB::PostProcessingImpl() {
  GetOutput() -= GetInput();
  return GetOutput() > 0;
}

}  // namespace dorofeev_i_bitwise_sort_double_eo_batcher_merge
