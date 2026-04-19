#pragma once
#include <vector>
#include "solonin_v_radix_sort_batcher/common/include/common.hpp"
#include "task/include/task.hpp"

namespace solonin_v_radix_sort_batcher {
class RadixSortBatcherSTL : public BaseTask {
 public:
  static constexpr ppc::task::TypeOfTask GetStaticTypeOfTask() {
    return ppc::task::TypeOfTask::kSTL;
  }
  explicit RadixSortBatcherSTL(const InType &in);

 private:
  bool ValidationImpl() override;
  bool PreProcessingImpl() override;
  bool RunImpl() override;
  bool PostProcessingImpl() override;

  static void SortByDigit(std::vector<int> &data, size_t pos);
  static void SortBlock(std::vector<int> &data, int lo, int hi);
  static std::vector<int> OddEvenMerge(const std::vector<int> &a, const std::vector<int> &b);
};
}  // namespace solonin_v_radix_sort_batcher
