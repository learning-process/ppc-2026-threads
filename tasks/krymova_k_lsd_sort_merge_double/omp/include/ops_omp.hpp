#pragma once

#include <cstdint>

#include "krymova_k_lsd_sort_merge_double/common/include/common.hpp"
#include "task/include/task.hpp"

namespace krymova_k_lsd_sort_merge_double {

class KrymovaKLsdSortMergeDoubleOMP : public BaseTask {
 public:
  static constexpr ppc::task::TypeOfTask GetStaticTypeOfTask() {
    return ppc::task::TypeOfTask::kOMP;
  }
  explicit KrymovaKLsdSortMergeDoubleOMP(const InType &in);

 private:
  bool ValidationImpl() override;
  bool PreProcessingImpl() override;
  bool RunImpl() override;
  bool PostProcessingImpl() override;

  static void LSDSortDouble(double *arr, int size);
  static void LSDSortDoubleParallel(double *arr, int size);
  static void IterativeMergeSort(double *arr, int size, int portion);
  static void MergeSections(double *left, const double *right, int left_size, int right_size);
  static void SortSectionsParallel(double *arr, int size, int portion);

  static uint64_t DoubleToULL(double d);
  static double ULLToDouble(uint64_t ull);

  int num_threads_;
};

}  // namespace krymova_k_lsd_sort_merge_double
