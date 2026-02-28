#pragma once

#include <vector>

#include "shekhirev_v_hoare_batcher_sort_seq/common/include/common.hpp"
#include "task/include/task.hpp"

namespace shekhirev_v_hoare_batcher_sort_seq {

class ShekhirevHoareBatcherSortSEQ : public BaseTask {
 public:
  static constexpr ppc::task::TypeOfTask GetStaticTypeOfTask() {
    return ppc::task::TypeOfTask::kSEQ;
  }
  explicit ShekhirevHoareBatcherSortSEQ(const InType &in);

 private:
  bool ValidationImpl() override;
  bool PreProcessingImpl() override;
  bool RunImpl() override;
  bool PostProcessingImpl() override;

  static void HoareSort(std::vector<int> &arr, int left, int right);
  static void BatcherMerge(std::vector<int> &arr, int left, int right, int r);
};

}  // namespace shekhirev_v_hoare_batcher_sort_seq
