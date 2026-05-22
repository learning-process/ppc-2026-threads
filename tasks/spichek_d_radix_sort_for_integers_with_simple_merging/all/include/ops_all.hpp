#pragma once

#include <vector>

#include "spichek_d_radix_sort_for_integers_with_simple_merging/common/include/common.hpp"
#include "task/include/task.hpp"

namespace spichek_d_radix_sort_for_integers_with_simple_merging {

class SpichekDRadixSortALL : public BaseTask {
 public:
  static constexpr ppc::task::TypeOfTask GetStaticTypeOfTask() {
    return ppc::task::TypeOfTask::kALL;
  }
  explicit SpichekDRadixSortALL(const InType &in);

 private:
  bool ValidationImpl() override;
  bool PreProcessingImpl() override;
  bool RunImpl() override;
  bool PostProcessingImpl() override;

  static void RadixSort(std::vector<int> &data);
  static std::vector<std::vector<int>> SplitData(const std::vector<int> &data, int num_parts);
  static std::vector<int> MergeData(std::vector<std::vector<int>> &parts);
};

}  // namespace spichek_d_radix_sort_for_integers_with_simple_merging
