#pragma once

#include <vector>

#include "task/include/task.hpp"
#include "yushkova_p_hoare_sorting_simple_merging/common/include/common.hpp"

namespace yushkova_p_hoare_sorting_simple_merging {

class YushkovaPHoareSortingSimpleMergingSTL : public BaseTask {
 public:
  static constexpr ppc::task::TypeOfTask GetStaticTypeOfTask() {
    return ppc::task::TypeOfTask::kSTL;
  }
  explicit YushkovaPHoareSortingSimpleMergingSTL(const InType &in);

 private:
  static int HoarePartition(std::vector<int> &values, int left, int right);
  static void HoareQuickSort(std::vector<int> &values, int left, int right);
  static std::vector<int> SimpleMerge(const std::vector<int> &left, const std::vector<int> &right);
  static void SortHalfIfNeeded(std::vector<int> &values);
  static void SortHalvesSequential(std::vector<int> &left, std::vector<int> &right);
  static void SortHalvesParallel(std::vector<int> &left, std::vector<int> &right);

  bool ValidationImpl() override;
  bool PreProcessingImpl() override;
  bool RunImpl() override;
  bool PostProcessingImpl() override;
};

}  // namespace yushkova_p_hoare_sorting_simple_merging
