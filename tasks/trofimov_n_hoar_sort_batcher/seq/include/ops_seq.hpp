#pragma once

#include "task/include/task.hpp"
#include "trofimov_n_hoar_sort_batcher/common/include/common.hpp"

namespace trofimov_n_hoar_sort_batcher {

class TrofimovNHoarSortBatcherSEQ : public BaseTask {
 public:
  static constexpr ppc::task::TypeOfTask GetStaticTypeOfTask() {
    return ppc::task::TypeOfTask::kSEQ;
  }
  explicit TrofimovNHoarSortBatcherSEQ(const InType &in);

 private:
  bool ValidationImpl() override;
  bool PreProcessingImpl() override;
  bool RunImpl() override;
  bool PostProcessingImpl() override;
};

}  // namespace trofimov_n_hoar_sort_batcher
