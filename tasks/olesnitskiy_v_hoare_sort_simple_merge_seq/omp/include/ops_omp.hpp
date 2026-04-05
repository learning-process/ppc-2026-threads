#pragma once

#include <vector>

#include "olesnitskiy_v_hoare_sort_simple_merge_seq/common/include/common.hpp"
#include "task/include/task.hpp"

namespace olesnitskiy_v_hoare_sort_simple_merge_seq {

class OlesnitskiyVHoareSortSimpleMergeOMP : public BaseTask {
 public:
  static constexpr ppc::task::TypeOfTask GetStaticTypeOfTask() {
    return ppc::task::TypeOfTask::kOMP;
  }
  explicit OlesnitskiyVHoareSortSimpleMergeOMP(const InType &in);

 private:
  static int HoarePartition(std::vector<int> &values, int left, int right);
  static void HoareQuickSort(std::vector<int> &values, int left, int right);
  static void Merge(std::vector<int> &values, int left, int mid, int right);

  bool ValidationImpl() override;
  bool PreProcessingImpl() override;
  bool RunImpl() override;
  bool PostProcessingImpl() override;
};

}  // namespace olesnitskiy_v_hoare_sort_simple_merge_seq
