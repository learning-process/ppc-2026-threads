#pragma once

#include "gonozov_l_bitwise_sorting_double_Batcher_merge/common/include/common.hpp"
#include "task/include/task.hpp"

namespace gonozov_l_bitwise_sorting_double_Batcher_merge {

class GonozovLBitSortBatcherMergeSEQ : public BaseTask {
 public:
  static constexpr ppc::task::TypeOfTask GetStaticTypeOfTask() {
    return ppc::task::TypeOfTask::kSEQ;
  }
  explicit GonozovLBitSortBatcherMergeSEQ(const InType &in);

 private:
  bool ValidationImpl() override;
  bool PreProcessingImpl() override;
  bool RunImpl() override;
  bool PostProcessingImpl() override;
};

}  // namespace gonozov_l_bitwise_sorting_double_Batcher_merge
