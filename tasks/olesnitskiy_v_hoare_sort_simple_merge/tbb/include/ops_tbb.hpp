#pragma once

#include <cstddef>
#include <vector>

#include "olesnitskiy_v_hoare_sort_simple_merge/common/include/common.hpp"
#include "task/include/task.hpp"

namespace olesnitskiy_v_hoare_sort_simple_merge {

class OlesnitskiyVHoareSortSimpleMergeTBB : public BaseTask {
 public:
  static constexpr ppc::task::TypeOfTask GetStaticTypeOfTask() {
    return ppc::task::TypeOfTask::kTBB;
  }
  explicit OlesnitskiyVHoareSortSimpleMergeTBB(const InType &in);

 private:
  static int HoarePartition(std::vector<int> &array, int left, int right);
  static void HoareQuickSort(std::vector<int> &array, int left, int right);
  static void SimpleMerge(const std::vector<int> &source, std::vector<int> &destination, size_t left, size_t middle,
                          size_t right);

  bool ValidationImpl() override;
  bool PreProcessingImpl() override;
  bool RunImpl() override;
  bool PostProcessingImpl() override;

  std::vector<int> data_;
};

}  // namespace olesnitskiy_v_hoare_sort_simple_merge
