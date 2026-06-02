#pragma once

#include <vector>

#include "krasnopevtseva_v_hoare_batcher_sort/common/include/common.hpp"
#include "task/include/task.hpp"

namespace krasnopevtseva_v_hoare_batcher_sort {

class KrasnopevtsevaVHoareBatcherSortSTL : public BaseTask {
 public:
  static constexpr ppc::task::TypeOfTask GetStaticTypeOfTask() {
    return ppc::task::TypeOfTask::kSTL;
  }
  explicit KrasnopevtsevaVHoareBatcherSortSTL(const InType &in);

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
  static void BatcherMergeBlocksStep(int *left_ptr, int &left_sz, int *right_ptr, int &right_sz);
  static void BatcherMergeLevel(int step, std::vector<Chunk> &chunks, int thread_input_size, int par_if_greater);
  static void ParallelSortChunks(std::vector<int> &arr, std::vector<Chunk> &chunks);
  static void SetupChunks(std::vector<Chunk> &chunks, int n, int numthreads);
  static int GetNumThreads(int n);
};
}  // namespace krasnopevtseva_v_hoare_batcher_sort
