#include "dorofeev_i_bitwise_sort_double_eo_batcher_merge/stl/include/ops_stl.hpp"

#include "dorofeev_i_bitwise_sort_double_eo_batcher_merge/common/include/common.hpp"

namespace dorofeev_i_bitwise_sort_double_eo_batcher_merge {

DorofeevIBitwiseSortDoubleEOBatcherMergeSTL::DorofeevIBitwiseSortDoubleEOBatcherMergeSTL(const InType &in) {
  SetTypeOfTask(GetStaticTypeOfTask());
  GetInput() = in;
}

bool DorofeevIBitwiseSortDoubleEOBatcherMergeSTL::ValidationImpl() {
  return true;
}
bool DorofeevIBitwiseSortDoubleEOBatcherMergeSTL::PreProcessingImpl() {
  return true;
}
bool DorofeevIBitwiseSortDoubleEOBatcherMergeSTL::RunImpl() {
  return true;
}
bool DorofeevIBitwiseSortDoubleEOBatcherMergeSTL::PostProcessingImpl() {
  return true;
}

}  // namespace dorofeev_i_bitwise_sort_double_eo_batcher_merge
