#pragma once

#include "dorofeev_i_bitwise_sort_double_eo_batcher_merge/common/include/common.hpp"
#include "task/include/task.hpp"

namespace dorofeev_i_bitwise_sort_double_eo_batcher_merge {

class DorofeevIBitwiseSortDoubleEOBatcherMergeSTL : public BaseTask {
 public:
  static constexpr ppc::task::TypeOfTask GetStaticTypeOfTask() {
    return ppc::task::TypeOfTask::kSTL;
  }
  explicit DorofeevIBitwiseSortDoubleEOBatcherMergeSTL(const InType &in);

 private:
  bool ValidationImpl() override;
  bool PreProcessingImpl() override;
  bool RunImpl() override;
  bool PostProcessingImpl() override;
};

}  // namespace dorofeev_i_bitwise_sort_double_eo_batcher_merge
