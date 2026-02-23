#pragma once

#include "krasnopevtseva_v_hoare_batcher_sort/common/include/common.hpp"
#include "task/include/task.hpp"

namespace krasnopevtseva_v_hoare_batcher_sort {

class KrasnopevtsevaVHoareBatcherSortSEQ : public BaseTask {
 public:
  static constexpr ppc::task::TypeOfTask GetStaticTypeOfTask() {
    return ppc::task::TypeOfTask::kSEQ;
  }
  explicit KrasnopevtsevaVHoareBatcherSortSEQ(const InType &in);

 private:
  bool ValidationImpl() override;
  bool PreProcessingImpl() override;
  bool RunImpl() override;
  void compareAndSwap(int &a, int &b);
  void batcherMerge(std::vector<int> &arr, int left, int right);
  void quickBatcherSort(std::vector<int> &arr, int left, int right);
  bool PostProcessingImpl() override;
};

}  // namespace krasnopevtseva_v_hoare_batcher_sort
