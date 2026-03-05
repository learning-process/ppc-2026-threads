#pragma once

#include "timofeev_n_radix_batcher_sort/common/include/common.hpp"
#include "task/include/task.hpp"

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

  int loggo(int inputa);
  void CompExch(int &a, int &b, int digit);
  void BubbleSort(std::vector<int> &arr, int digit, int left, int right);
  void ComparR(int &a, int &b);
  void OddEvenMerge(std::vector<int> &arr, size_t lft, size_t n);

};

}  // namespace timofeev_n_radix_batcher_sort_threads
