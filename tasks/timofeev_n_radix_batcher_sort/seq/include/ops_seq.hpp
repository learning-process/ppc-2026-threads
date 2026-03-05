#pragma once

#include <vector>

#include "task/include/task.hpp"
#include "timofeev_n_radix_batcher_sort/common/include/common.hpp"

namespace timofeev_n_radix_batcher_sort_threads {

class TimofeevNRadixBatcherSEQ : public BaseTask {
 public:
  static constexpr ppc::task::TypeOfTask GetStaticTypeOfTask() {
    return ppc::task::TypeOfTask::kSEQ;
  }
  explicit TimofeevNRadixBatcherSEQ(const InType &in);

 private:
  bool ValidationImpl() override;
  bool PreProcessingImpl() override;
  bool RunImpl() override;
  bool PostProcessingImpl() override;

  static int Loggo(int inputa);
  static void CompExch(int &a, int &b, int digit);
  static void BubbleSort(std::vector<int> &arr, int digit, int left, int right);
  static void ComparR(int &a, int &b);
  static void OddEvenMerge(std::vector<int> &arr, int lft, int n);
};

}  // namespace timofeev_n_radix_batcher_sort_threads
