#pragma once

#include <vector>

#include "nikitina_v_hoar_sort_batcher/common/include/common.hpp"
#include "task/include/task.hpp"

namespace nikitina_v_hoar_sort_batcher {

class HoareSortBatcherSEQ : public BaseTask {
 public:
  static constexpr ppc::task::TypeOfTask GetStaticTypeOfTask() {
    return ppc::task::TypeOfTask::kSEQ;
  }

  explicit HoareSortBatcherSEQ(const InType &in);

 private:
  bool ValidationImpl() override;
  bool PreProcessingImpl() override;
  bool RunImpl() override;
  bool PostProcessingImpl() override;
};

}  // namespace nikitina_v_hoar_sort_batcher
