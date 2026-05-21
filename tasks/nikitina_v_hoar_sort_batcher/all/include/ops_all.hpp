#pragma once

#include "nikitina_v_hoar_sort_batcher/common/include/common.hpp"
#include "task/include/task.hpp"

namespace nikitina_v_hoar_sort_batcher {

class HoareSortBatcherALL : public BaseTask {
 public:
  static constexpr ppc::task::TypeOfTask GetStaticTypeOfTask() {
    return ppc::task::TypeOfTask::kALL;
  }
  explicit HoareSortBatcherALL(const InType &in);

 private:
  bool ValidationImpl() override;
  bool PreProcessingImpl() override;
  bool RunImpl() override;
  bool PostProcessingImpl() override;

  OutType data_;
};

}  // namespace nikitina_v_hoar_sort_batcher
