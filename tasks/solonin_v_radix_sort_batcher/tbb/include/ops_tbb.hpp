#pragma once
#include <cstddef>
#include <vector>

#include "solonin_v_radix_sort_batcher/common/include/common.hpp"
#include "task/include/task.hpp"

namespace solonin_v_radix_sort_batcher {
class RadixSortBatcherTBB : public BaseTask {
 public:
  static constexpr ppc::task::TypeOfTask GetStaticTypeOfTask() {
    return ppc::task::TypeOfTask::kTBB;
  }
  explicit RadixSortBatcherTBB(const InType &in);

 private:
  bool ValidationImpl() override;
  bool PreProcessingImpl() override;
  bool RunImpl() override;
  bool PostProcessingImpl() override;

  static void SortByDigit(std::vector<int> &data, size_t pos);
  static void SortSegment(std::vector<int> &data, int lo, int hi);
  static std::vector<int> MergeBatcher(const std::vector<int> &left, const std::vector<int> &right);
};
}  // namespace solonin_v_radix_sort_batcher
