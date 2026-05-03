#pragma once

#include "akimov_i_radixsort_int_merge/common/include/common.hpp"
#include "task/include/task.hpp"

namespace akimov_i_radixsort_int_merge {

class AkimovIRadixSortIntMergeOMP : public BaseTask {
 public:
  static constexpr ppc::task::TypeOfTask GetStaticTypeOfTask() {
    return ppc::task::TypeOfTask::kOMP;
  }
  explicit AkimovIRadixSortIntMergeOMP(const InType &in);

 private:
  bool ValidationImpl() override;
  bool PreProcessingImpl() override;
  bool RunImpl() override;
  bool PostProcessingImpl() override;
};

}  // namespace akimov_i_radixsort_int_merge
