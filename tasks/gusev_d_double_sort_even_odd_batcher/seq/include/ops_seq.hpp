#pragma once

#include <vector>

#include "task/include/task.hpp"

namespace gusev_d_double_sort_even_odd_batcher_task_threads {

using InType = std::vector<double>;
using OutType = std::vector<double>;
using BaseTask = ppc::task::Task<InType, OutType>;

class DoubleSortEvenOddBatcherSEQ : public BaseTask {
 public:
  static constexpr ppc::task::TypeOfTask GetStaticTypeOfTask() {
    return ppc::task::TypeOfTask::kSEQ;
  }
  explicit DoubleSortEvenOddBatcherSEQ(const InType &in);

 private:
  bool ValidationImpl() override;
  bool PreProcessingImpl() override;
  bool RunImpl() override;
  bool PostProcessingImpl() override;

  InType input_data_;
  OutType result_data_;
};

}  // namespace gusev_d_double_sort_even_odd_batcher_task_threads
