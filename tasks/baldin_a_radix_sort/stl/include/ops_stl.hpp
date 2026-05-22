#pragma once

#include "baldin_a_radix_sort/common/include/common.hpp"
#include "task/include/task.hpp"

namespace baldin_a_radix_sort {

class BaldinARadixSortSTL : public BaseTask {
 public:
  static constexpr ppc::task::TypeOfTask GetStaticTypeOfTask() {
    return ppc::task::TypeOfTask::kSTL;
  }
  explicit BaldinARadixSortSTL(const InType &in);

 private:
  bool ValidationImpl() override;
  bool PreProcessingImpl() override;
  bool RunImpl() override;
  bool PostProcessingImpl() override;
};

}  // namespace baldin_a_radix_sort
