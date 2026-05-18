#pragma once

#include <vector>

#include "olesnitskiy_v_hoare_sort_simple_merge/common/include/common.hpp"
#include "task/include/task.hpp"

namespace olesnitskiy_v_hoare_sort_simple_merge {

class OlesnitskiyVHoareSortSimpleMergeOMP : public BaseTask {
 public:
  static constexpr ppc::task::TypeOfTask GetStaticTypeOfTask() {
    return ppc::task::TypeOfTask::kOMP;
  }
  explicit OlesnitskiyVHoareSortSimpleMergeOMP(const InType &in);

 private:
  static int HoarePartition(std::vector<int> &array, int left, int right);
  static void HoareQuickSort(std::vector<int> &array, int left, int right);
  static std::vector<int> SimpleMerge(const std::vector<int> &left_part, const std::vector<int> &right_part);

  bool ValidationImpl() override;
  bool PreProcessingImpl() override;
  bool RunImpl() override;
  bool PostProcessingImpl() override;

  std::vector<int> data_;
};

}  // namespace olesnitskiy_v_hoare_sort_simple_merge
