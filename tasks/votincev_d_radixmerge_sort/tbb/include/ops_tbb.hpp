#pragma once

#include <cstdint>
#include <vector>

#include "task/include/task.hpp"
#include "votincev_d_radixmerge_sort/common/include/common.hpp"

namespace votincev_d_radixmerge_sort {

class VotincevDRadixMergeSortTBB : public BaseTask {
 public:
  static constexpr ppc::task::TypeOfTask GetStaticTypeOfTask() {
    return ppc::task::TypeOfTask::kTBB;
  }
  explicit VotincevDRadixMergeSortTBB(const InType &in);

 private:
  bool ValidationImpl() override;
  bool PreProcessingImpl() override;
  bool RunImpl() override;
  bool PostProcessingImpl() override;

  // ==============================
  // мои дополнительные функции ===
  static void LocalRadixSort(uint32_t *begin, uint32_t *end);
  static void Merge(uint32_t *data, int32_t left, int32_t mid, int32_t right, uint32_t *temp);
  static void ParallelRadixMergeSort(uint32_t *data, int32_t left, int32_t right, uint32_t *temp);
};
}  // namespace votincev_d_radixmerge_sort
