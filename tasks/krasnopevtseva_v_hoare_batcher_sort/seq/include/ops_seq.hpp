#pragma once

#include <vector>

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
  static void CompareAndSwap(int &a, int &b);
  static void BatcherMerge(std::vector<int> &arr, int left, int right);
  static void QuickBatcherSort(std::vector<int> &arr, int left, int right);
  static void InsertionSort(std::vector<int> &arr, int left, int right);
  static int Partition(std::vector<int> &arr, int left, int right);
  bool PostProcessingImpl() override;
};

}  // namespace krasnopevtseva_v_hoare_batcher_sort
