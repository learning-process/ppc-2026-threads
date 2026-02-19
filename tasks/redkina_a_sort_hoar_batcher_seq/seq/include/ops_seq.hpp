// ops_seq.hpp
#pragma once

#include "redkina_a_sort_hoar_batcher_seq/common/include/common.hpp"
#include "task/include/task.hpp"

namespace redkina_a_sort_hoar_batcher_seq {

class RedkinaASortHoarBatcherSEQ : public BaseTask {
 public:
  static constexpr ppc::task::TypeOfTask GetStaticTypeOfTask() {
    return ppc::task::TypeOfTask::kSEQ;
  }
  explicit RedkinaASortHoarBatcherSEQ(const InType& in);

 private:
  bool ValidationImpl() override;
  bool PreProcessingImpl() override;
  bool RunImpl() override;
  bool PostProcessingImpl() override;
};

}  // namespace redkina_a_sort_hoar_batcher_seq