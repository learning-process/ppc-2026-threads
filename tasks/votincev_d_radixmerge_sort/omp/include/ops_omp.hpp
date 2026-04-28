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
};
}  // namespace votincev_d_radixmerge_sort
