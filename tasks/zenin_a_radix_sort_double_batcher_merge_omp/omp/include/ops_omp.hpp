#pragma once

#include "task/include/task.hpp"
#include "zenin_a_radix_sort_double_batcher_merge_omp/common/include/common.hpp"

namespace zenin_a_radix_sort_double_batcher_merge_omp {

class ZeninARadixSortDoubleBatcherMergeOMP : public BaseTask {
 public:
  static constexpr ppc::task::TypeOfTask GetStaticTypeOfTask() {
    return ppc::task::TypeOfTask::kOMP;
  }
  explicit ZeninARadixSortDoubleBatcherMergeOMP(const InType &in);

 private:
  bool ValidationImpl() override;
  bool PreProcessingImpl() override;
  bool RunImpl() override;
  bool PostProcessingImpl() override;
};

}  // namespace zenin_a_radix_sort_double_batcher_merge_omp
