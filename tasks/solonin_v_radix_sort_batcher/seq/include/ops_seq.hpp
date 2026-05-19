#pragma once
#include "solonin_v_radix_sort_batcher/common/include/common.hpp"
#include "task/include/task.hpp"
#include <cstddef>
#include <vector>

namespace solonin_v_radix_sort_batcher {
class RadixSortBatcherSEQ : public BaseTask {
 public:
  static constexpr ppc::task::TypeOfTask GetStaticTypeOfTask() {
    return ppc::task::TypeOfTask::kSEQ;
  }
  explicit RadixSortBatcherSEQ(const InType &in);

 private:
  bool ValidationImpl() override;
  bool PreProcessingImpl() override;
  bool RunImpl() override;
  bool PostProcessingImpl() override;

  static void SortByDigit(std::vector<int> &data, size_t pos);
  static void RadixSort(std::vector<int> &data);
};
}  // namespace solonin_v_radix_sort_batcher
