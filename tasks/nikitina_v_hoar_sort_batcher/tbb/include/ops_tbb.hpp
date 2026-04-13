#pragma once

#include "nikitina_v_hoar_sort_batcher/common/include/common.hpp"
#include "task/include/task.hpp"
namespace nikitina_v_hoar_sort_batcher {

class HoareSortBatcherTBB : public BaseTask {
 public:
  static constexpr ppc::task::TypeOfTask GetStaticTypeOfTask() {
    return ppc::task::TypeOfTask::kTBB;
  }

  explicit HoareSortBatcherTBB(const InType &in);

 private:
  bool ValidationImpl() override;
  bool PreProcessingImpl() override;
  bool RunImpl() override;
  bool PostProcessingImpl() override;
};

}  // namespace nikitina_v_hoar_sort_batcher
