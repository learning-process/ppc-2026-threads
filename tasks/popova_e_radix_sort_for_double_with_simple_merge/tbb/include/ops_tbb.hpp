#pragma once

#include <vector>

#include "popova_e_radix_sort_for_double_with_simple_merge/common/include/common.hpp"
#include "task/include/task.hpp"

namespace popova_e_radix_sort_for_double_with_simple_merge_threads {

class PopovaERadixSorForDoubleWithSimpleMergeTBB : public BaseTask {
 public:
  static constexpr ppc::task::TypeOfTask GetStaticTypeOfTask() {
    return ppc::task::TypeOfTask::kTBB;
  }
  explicit PopovaERadixSorForDoubleWithSimpleMergeTBB(const InType &in);

 private:
  bool ValidationImpl() override;
  bool PreProcessingImpl() override;
  bool RunImpl() override;
  bool PostProcessingImpl() override;

  std::vector<double> array_;   // массив для сортировки
  std::vector<double> result_;  // результат сортировки
};

}  // namespace popova_e_radix_sort_for_double_with_simple_merge_threads
