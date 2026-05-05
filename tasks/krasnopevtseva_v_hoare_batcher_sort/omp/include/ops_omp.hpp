#pragma once

#include <vector>

#include "krasnopevtseva_v_hoare_batcher_sort/common/include/common.hpp"
#include "task/include/task.hpp"

namespace krasnopevtseva_v_hoare_batcher_sort {

class KrasnopevtsevaVHoareBatcherSortOMP : public BaseTask {
 public:
  static constexpr ppc::task::TypeOfTask GetStaticTypeOfTask() {
    return ppc::task::TypeOfTask::kOMP;
  }
  explicit KrasnopevtsevaVHoareBatcherSortOMP(const InType &in);

 private:
  bool ValidationImpl() override;
  bool PreProcessingImpl() override;
  bool RunImpl() override;
  bool PostProcessingImpl() override;
  static int Partition(std::vector<int> &arr, int first, int last);
  static void QuickSort(std::vector<int> &arr, int first, int last);
  static void InsertionSort(std::vector<int> &arr, int first, int last);
  static void BatcherMergeBlocksStep(int *left_pointer, int &left_size, int *right_pointer, int &right_size);
  static void BatcherMerge(int thread_input_size, std::vector<int *> &pointers, std::vector<int> &sizes,
                           int par_if_greater);
};
}  // namespace krasnopevtseva_v_hoare_batcher_sort
