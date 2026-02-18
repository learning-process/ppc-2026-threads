#pragma once

#include "zenin_a_radix_sort_double_batcher_merge_seq/common/include/common.hpp"
#include "task/include/task.hpp"

namespace zenin_a_radix_sort_double_batcher_merge_seq {

class ZeninARadixSortDoubleBatcherMerge_SEQSEQ : public BaseTask {
 public:
  static constexpr ppc::task::TypeOfTask GetStaticTypeOfTask() {
    return ppc::task::TypeOfTask::kSEQ;
  }
  explicit ZeninARadixSortDoubleBatcherMerge_SEQSEQ(const InType &in);

 private:
  bool ValidationImpl() override;
  bool PreProcessingImpl() override;
  bool RunImpl() override;
  bool PostProcessingImpl() override;
};

}  // namespace zenin_a_radix_sort_double_batcher_merge_seq
