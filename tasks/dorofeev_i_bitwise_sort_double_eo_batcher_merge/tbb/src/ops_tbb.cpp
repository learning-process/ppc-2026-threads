#include "dorofeev_i_bitwise_sort_double_eo_batcher_merge/tbb/include/ops_tbb.hpp"

#include "dorofeev_i_bitwise_sort_double_eo_batcher_merge/common/include/common.hpp"

namespace dorofeev_i_bitwise_sort_double_eo_batcher_merge {

DorofeevIBitwiseSortDoubleEOBatcherMergeTBB::DorofeevIBitwiseSortDoubleEOBatcherMergeTBB(const InType &in) {
  SetTypeOfTask(GetStaticTypeOfTask());
  GetInput() = in;
}

bool DorofeevIBitwiseSortDoubleEOBatcherMergeTBB::ValidationImpl() {
  return true;
}
bool DorofeevIBitwiseSortDoubleEOBatcherMergeTBB::PreProcessingImpl() {
  return true;
}
bool DorofeevIBitwiseSortDoubleEOBatcherMergeTBB::RunImpl() {
  return true;
}
bool DorofeevIBitwiseSortDoubleEOBatcherMergeTBB::PostProcessingImpl() {
  return true;
}

}  // namespace dorofeev_i_bitwise_sort_double_eo_batcher_merge
