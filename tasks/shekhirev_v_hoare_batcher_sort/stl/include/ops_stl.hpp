#pragma once

#include "shekhirev_v_hoare_batcher_sort/common/include/common.hpp"
#include "task/include/task.hpp"

namespace shekhirev_v_hoare_batcher_sort {

class ShekhirevHoareBatcherSortSTL : public BaseTask {
 public:
  static constexpr ppc::task::TypeOfTask GetStaticTypeOfTask() {
    return ppc::task::TypeOfTask::kSTL;
  }

  explicit ShekhirevHoareBatcherSortSTL(const InType &in);

 private:
  bool ValidationImpl() override;
  bool PreProcessingImpl() override;
  bool RunImpl() override;
  bool PostProcessingImpl() override;
};

}  // namespace shekhirev_v_hoare_batcher_sort
