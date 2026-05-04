#pragma once

#include <stack>
#include <utility>

#include "rozenberg_a_quicksort_simple_merge/common/include/common.hpp"
#include "task/include/task.hpp"

namespace rozenberg_a_quicksort_simple_merge {

class RozenbergAQuicksortSimpleMergeTBB : public BaseTask {
 public:
  static constexpr ppc::task::TypeOfTask GetStaticTypeOfTask() {
    return ppc::task::TypeOfTask::kTBB;
  }
  explicit RozenbergAQuicksortSimpleMergeTBB(const InType &in);

 private:
  bool ValidationImpl() override;
  bool PreProcessingImpl() override;
  bool RunImpl() override;
  bool PostProcessingImpl() override;
  static std::pair<int, int> Partition(InType &data, int left, int right);
  static void PushSubarrays(std::stack<std::pair<int, int>> &stack, int left, int right, int i, int j);
  static void Quicksort(InType &data, int low, int high);
  static void Merge(InType &data, int left, int mid, int right);
};

}  // namespace rozenberg_a_quicksort_simple_merge
