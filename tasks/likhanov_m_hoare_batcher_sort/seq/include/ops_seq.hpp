#pragma once

#include "likhanov_m_hoare_batcher_sort/common/include/common.hpp"
#include "task/include/task.hpp"

namespace likhanov_m_hoare_batcher_sort {

class LikhanovMHoareBatcherSortSEQ : public BaseTask {
 public:
  static constexpr ppc::task::TypeOfTask GetStaticTypeOfTask() {
    return ppc::task::TypeOfTask::kSEQ;
  }
  explicit LikhanovMHoareBatcherSortSEQ(const InType &in);

 private:
  bool ValidationImpl() override;
  bool PreProcessingImpl() override;
  bool RunImpl() override;
  bool PostProcessingImpl() override;
};

}  // namespace likhanov_m_hoare_batcher_sort
