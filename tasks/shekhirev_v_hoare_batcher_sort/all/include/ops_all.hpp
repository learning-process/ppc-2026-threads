#pragma once

#include "shekhirev_v_hoare_batcher_sort/common/include/common.hpp"
#include "task/include/task.hpp"

namespace shekhirev_v_hoare_batcher_sort {

class ShekhirevHoareBatcherSortALL : public BaseTask {
 public:
  static constexpr ppc::task::TypeOfTask GetStaticTypeOfTask() {
    return ppc::task::TypeOfTask::kALL;
  }

  explicit ShekhirevHoareBatcherSortALL(const InType &in);

 private:
  bool ValidationImpl() override;
  bool PreProcessingImpl() override;
  bool RunImpl() override;
  bool PostProcessingImpl() override;

  InType input_;
  OutType output_;
};

}  // namespace shekhirev_v_hoare_batcher_sort
