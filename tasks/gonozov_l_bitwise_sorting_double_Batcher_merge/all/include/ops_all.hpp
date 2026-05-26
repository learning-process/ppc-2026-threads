#pragma once

#include <vector>

#include "gonozov_l_bitwise_sorting_double_Batcher_merge/common/include/common.hpp"
#include "task/include/task.hpp"

namespace gonozov_l_bitwise_sorting_double_batcher_merge {

class GonozovLBitSortBatcherMergeALL : public BaseTask {
 public:
  static constexpr ppc::task::TypeOfTask GetStaticTypeOfTask() {
    return ppc::task::TypeOfTask::kALL;
  }
  explicit GonozovLBitSortBatcherMergeALL(const InType &in);

 private:
  bool ValidationImpl() override;
  bool PreProcessingImpl() override;
  bool RunImpl() override;
  bool PostProcessingImpl() override;

  std::vector<double> local_data_;
};

}  // namespace gonozov_l_bitwise_sorting_double_batcher_merge
