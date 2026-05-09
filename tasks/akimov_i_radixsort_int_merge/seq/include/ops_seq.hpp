#pragma once

#include "akimov_i_radixsort_int_merge/common/include/common.hpp"
#include "task/include/task.hpp"

namespace akimov_i_radixsort_int_merge {

class AkimovIRadixSortIntMergeSEQ : public BaseTask {
 public:
  static constexpr ppc::task::TypeOfTask GetStaticTypeOfTask() {
    return ppc::task::TypeOfTask::kSEQ;
  }
  explicit AkimovIRadixSortIntMergeSEQ(const InType &in);

 private:
  bool ValidationImpl() override;
  bool PreProcessingImpl() override;
  bool RunImpl() override;
  bool PostProcessingImpl() override;
};

}  // namespace akimov_i_radixsort_int_merge
