#pragma once

#include "popova_e_radix_sort_for_double_with_simple_merge/common/include/common.hpp"
#include "task/include/task.hpp"

namespace popova_e_radix_sort_for_double_with_simple_merge_threads {

class PopovaERadixSorForDoubleWithSimpleMergeSEQ : public BaseTask {
 public:
  static constexpr ppc::task::TypeOfTask GetStaticTypeOfTask() {
    return ppc::task::TypeOfTask::kSEQ;
  }
  explicit PopovaERadixSorForDoubleWithSimpleMergeSEQ(const InType &in);

 private:
  bool ValidationImpl() override;
  bool PreProcessingImpl() override;
  bool RunImpl() override;
  bool PostProcessingImpl() override;
};

}  // namespace popova_e_radix_sort_for_double_with_simple_merge_threads