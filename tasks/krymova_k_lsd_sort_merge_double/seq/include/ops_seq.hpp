#pragma once

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
  
  void LSDSortDouble(double* arr, int size);
  void IterativeMergeSort(double* arr, double* tmp, int size, int portion);
  
  static unsigned long long DoubleToULL(double d);
  static double ULLToDouble(unsigned long long ull);
};

}  // namespace krymova_k_lsd_sort_merge_double
