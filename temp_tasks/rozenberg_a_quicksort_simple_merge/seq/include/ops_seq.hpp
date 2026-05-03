#pragma once

#include <stack>
#include <utility>

#include "rozenberg_a_quicksort_simple_merge/common/include/common.hpp"
#include "task/include/task.hpp"

namespace rozenberg_a_quicksort_simple_merge {

class RozenbergAQuicksortSimpleMergeSEQ : public BaseTask {
 public:
  static constexpr ppc::task::TypeOfTask GetStaticTypeOfTask() {
    return ppc::task::TypeOfTask::kSEQ;
  }
  explicit RozenbergAQuicksortSimpleMergeSEQ(const InType &in);

 private:
  bool ValidationImpl() override;
  bool PreProcessingImpl() override;
  bool RunImpl() override;
  bool PostProcessingImpl() override;
  static std::pair<int, int> Partition(InType &data, int left, int right);
  static void PushSubarrays(std::stack<std::pair<int, int>> &stack, int left, int right, int i, int j);
  static void Quicksort(InType &data);
};

}  // namespace rozenberg_a_quicksort_simple_merge
