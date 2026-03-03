#include "zenin_a_radix_sort_double_batcher_merge_omp/omp/include/ops_omp.hpp"

#include <atomic>
#include <numeric>
#include <vector>

#include "util/include/util.hpp"
#include "zenin_a_radix_sort_double_batcher_merge_omp/common/include/common.hpp"

namespace zenin_a_radix_sort_double_batcher_merge_omp {

ZeninARadixSortDoubleBatcherMergeOMP::ZeninARadixSortDoubleBatcherMergeOMP(const InType &in) {
  SetTypeOfTask(GetStaticTypeOfTask());
  GetInput() = in;
  GetOutput() = {};
}

bool ZeninARadixSortDoubleBatcherMergeOMP::ValidationImpl() {
  return true;
}

bool ZeninARadixSortDoubleBatcherMergeOMP::PreProcessingImpl() {
  return true;
}

bool ZeninARadixSortDoubleBatcherMergeOMP::RunImpl() {
  // #pragma omp parallel default(none) shared(counter) num_threads(ppc::util::GetNumThreads())

  return true;
}

bool ZeninARadixSortDoubleBatcherMergeOMP::PostProcessingImpl() {
  return true;
}

}  // namespace zenin_a_radix_sort_double_batcher_merge_omp
