#pragma once

#include <cstdint>

#include "krymova_k_lsd_sort_merge_double/common/include/common.hpp"
#include "task/include/task.hpp"

namespace krymova_k_lsd_sort_merge_double {

class KrymovaKLsdSortMergeDoubleSTL : public BaseTask {
 public:
  static constexpr ppc::task::TypeOfTask GetStaticTypeOfTask() {
    return ppc::task::TypeOfTask::kSTL;
  }
  explicit KrymovaKLsdSortMergeDoubleSTL(const InType &in);

 private:
  bool ValidationImpl() override;
  bool PreProcessingImpl() override;
  bool RunImpl() override;
  bool PostProcessingImpl() override;

  static void LSDSortDoubleSequential(double *arr, int size);
  static void IterativeMergeSort(double *arr, int size, int portion, int num_threads);
  static void MergeSections(double *left, const double *right, int left_size, int right_size);
  static void SortSectionsParallel(double *arr, int size, int portion, int num_threads);

  static uint64_t DoubleToULL(double d);
  static double ULLToDouble(uint64_t ull);
};

}  // namespace krymova_k_lsd_sort_merge_double
