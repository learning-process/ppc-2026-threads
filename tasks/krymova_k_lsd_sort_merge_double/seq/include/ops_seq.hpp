#pragma once

#include <cstdint>

#include "krymova_k_lsd_sort_merge_double/common/include/common.hpp"
#include "task/include/task.hpp"

namespace krymova_k_lsd_sort_merge_double {

class KrymovaKLsdSortMergeDoubleSEQ : public BaseTask {
 public:
  static constexpr ppc::task::TypeOfTask GetStaticTypeOfTask() {
    return ppc::task::TypeOfTask::kSEQ;
  }
  explicit KrymovaKLsdSortMergeDoubleSEQ(const InType &in);

 private:
  bool ValidationImpl() override;
  bool PreProcessingImpl() override;
  bool RunImpl() override;
  bool PostProcessingImpl() override;

  static void LSDSortDouble(double* arr, int size);
  static void IterativeMergeSort(double* arr, double* tmp, int size, int portion);

  static uint64_t DoubleToULL(double d);
  static double ULLToDouble(uint64_t ull);
};

}  // namespace krymova_k_lsd_sort_merge_double