#pragma once

#include <cstdint>
#include <vector>

#include "task/include/task.hpp"
#include "votincev_d_radixmerge_sort/common/include/common.hpp"

namespace votincev_d_radixmerge_sort {

class VotincevDRadixMergeSortOMP : public BaseTask {
 public:
  static constexpr ppc::task::TypeOfTask GetStaticTypeOfTask() {
    return ppc::task::TypeOfTask::kOMP;
  }
  explicit VotincevDRadixMergeSortOMP(const InType &in);

 private:
  bool ValidationImpl() override;
  bool PreProcessingImpl() override;
  bool RunImpl() override;
  bool PostProcessingImpl() override;

  // ==============================
  // мои дополнительные функции ===
  static void SortByDigit(std::vector<int32_t> &array, int32_t exp);
  static void LocalRadixSort(uint32_t *begin, uint32_t *end);
  static void Merge(uint32_t *src, uint32_t *dst, int32_t left, int32_t mid, int32_t right);
};
}  // namespace votincev_d_radixmerge_sort
