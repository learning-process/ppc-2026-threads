#pragma once

#include "gonozov_l_bitwise_sorting_double_Batcher_merge/common/include/common.hpp"
#include "task/include/task.hpp"

namespace gonozov_l_bitwise_sorting_double_batcher_merge {

class GonozovLBitSortBatcherMergeTBB : public BaseTask {
 public:
  static constexpr ppc::task::TypeOfTask GetStaticTypeOfTask() {
    return ppc::task::TypeOfTask::kTBB;
  }
  explicit GonozovLBitSortBatcherMergeTBB(const InType &in);

 private:
  bool ValidationImpl() override;
  bool PreProcessingImpl() override;
  bool RunImpl() override;
  bool PostProcessingImpl() override;
};

}  // namespace gonozov_l_bitwise_sorting_double_batcher_merge
