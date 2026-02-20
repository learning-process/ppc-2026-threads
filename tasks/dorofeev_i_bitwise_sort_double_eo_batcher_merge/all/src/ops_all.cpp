#include "dorofeev_i_bitwise_sort_double_eo_batcher_merge/all/include/ops_all.hpp"

namespace dorofeev_i_bitwise_sort_double_eo_batcher_merge {

DorofeevIBitwiseSortDoubleEOBatcherMergeALL::DorofeevIBitwiseSortDoubleEOBatcherMergeALL(const InType &in) {
  SetTypeOfTask(GetStaticTypeOfTask());
  GetInput() = in;
}

bool DorofeevIBitwiseSortDoubleEOBatcherMergeALL::ValidationImpl() { return true; }
bool DorofeevIBitwiseSortDoubleEOBatcherMergeALL::PreProcessingImpl() { return true; }
bool DorofeevIBitwiseSortDoubleEOBatcherMergeALL::RunImpl() { return true; }
bool DorofeevIBitwiseSortDoubleEOBatcherMergeALL::PostProcessingImpl() { return true; }

}  // namespace dorofeev_i_bitwise_sort_double_eo_batcher_merge