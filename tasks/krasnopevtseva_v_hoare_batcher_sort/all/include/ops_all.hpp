#pragma once

#include <vector>

#include "krasnopevtseva_v_hoare_batcher_sort/common/include/common.hpp"
#include "task/include/task.hpp"

namespace krasnopevtseva_v_hoare_batcher_sort {

class KrasnopevtsevaVHoareBatcherSortALL : public BaseTask {
 public:
  static constexpr ppc::task::TypeOfTask GetStaticTypeOfTask() {
    return ppc::task::TypeOfTask::kALL;
  }
  explicit KrasnopevtsevaVHoareBatcherSortALL(const InType &in);

 private:
  bool ValidationImpl() override;
  bool PreProcessingImpl() override;
  bool RunImpl() override;
  bool PostProcessingImpl() override;

  struct Chunk {
    int *ptr;
    int size;
    int left;
    int right;
  };

  static int Partition(std::vector<int> &arr, int first, int last);
  static void InsertionSort(std::vector<int> &arr, int first, int last);
  static void QuickSort(std::vector<int> &arr, int first, int last);
  static void BatcherMergeBlocksStep(int *left_pointer, int &left_size, int *right_pointer, int &right_size);
  static void BatcherMerge(int thread_input_size, std::vector<int *> &pointers, std::vector<int> &sizes,
                           int par_if_greater);
  static void ParallelSortChunksOpenMP(std::vector<int> &res, int n, int numthreads);
  static void ParallelSortChunks(std::vector<int> &arr, int n, int numthreads);
  static void SortLocalData(std::vector<int> &data);
};
}  // namespace krasnopevtseva_v_hoare_batcher_sort
