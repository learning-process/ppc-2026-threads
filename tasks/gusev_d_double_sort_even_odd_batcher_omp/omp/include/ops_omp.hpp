#pragma once

#include "gusev_d_double_sort_even_odd_batcher_omp/common/include/common.hpp"
#include "task/include/task.hpp"

namespace gusev_d_double_sort_even_odd_batcher_omp_task_threads {

class DoubleSortEvenOddBatcherOMP : public BaseTask {
 public:
  static constexpr ppc::task::TypeOfTask GetStaticTypeOfTask() {
    return ppc::task::TypeOfTask::kOMP;
  }

  explicit DoubleSortEvenOddBatcherOMP(const InType &in);

 private:
  bool ValidationImpl() override;
  bool PreProcessingImpl() override;
  bool RunImpl() override;
  bool PostProcessingImpl() override;

  InType input_data_;
  OutType result_data_;
};

}  // namespace gusev_d_double_sort_even_odd_batcher_omp_task_threads
