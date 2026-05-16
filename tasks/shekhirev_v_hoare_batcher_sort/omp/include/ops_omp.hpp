#pragma once

#include "shekhirev_v_hoare_batcher_sort/common/include/common.hpp"
#include "task/include/task.hpp"

namespace shekhirev_v_hoare_batcher_sort {

class ShekhirevHoareBatcherSortOMP : public BaseTask {
 public:
  static constexpr ppc::task::TypeOfTask GetStaticTypeOfTask() {
    return ppc::task::TypeOfTask::kOMP;
  }

  explicit ShekhirevHoareBatcherSortOMP(const InType &in);

 protected:
  bool ValidationImpl() override;
  bool PreProcessingImpl() override;
  bool RunImpl() override;
  bool PostProcessingImpl() override;

 private:
  InType input_;
  OutType res_;
};

}  // namespace shekhirev_v_hoare_batcher_sort
