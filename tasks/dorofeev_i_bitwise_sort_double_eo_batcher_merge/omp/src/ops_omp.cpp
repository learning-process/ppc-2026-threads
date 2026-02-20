#include "dorofeev_i_bitwise_sort_double_eo_batcher_merge/omp/include/ops_omp.hpp"

#include "dorofeev_i_bitwise_sort_double_eo_batcher_merge/common/include/common.hpp"

namespace dorofeev_i_bitwise_sort_double_eo_batcher_merge {

DorofeevIBitwiseSortDoubleEOBatcherMergeOMP::DorofeevIBitwiseSortDoubleEOBatcherMergeOMP(const InType &in) {
  SetTypeOfTask(GetStaticTypeOfTask());
  GetInput() = in;
}

bool DorofeevIBitwiseSortDoubleEOBatcherMergeOMP::ValidationImpl() {
  return true;
}
bool DorofeevIBitwiseSortDoubleEOBatcherMergeOMP::PreProcessingImpl() {
  return true;
}
bool DorofeevIBitwiseSortDoubleEOBatcherMergeOMP::RunImpl() {
  return true;
}
bool DorofeevIBitwiseSortDoubleEOBatcherMergeOMP::PostProcessingImpl() {
  return true;
}

}  // namespace dorofeev_i_bitwise_sort_double_eo_batcher_merge
