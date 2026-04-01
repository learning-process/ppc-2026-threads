#pragma once

#include <cstddef>
#include <vector>

#include "task/include/task.hpp"
#include "timofeev_n_radix_batcher_sort/common/include/common.hpp"

namespace timofeev_n_radix_batcher_sort_threads {

class TimofeevNRadixBatcherTBB : public BaseTask {
 public:
  static constexpr ppc::task::TypeOfTask GetStaticTypeOfTask() {
    return ppc::task::TypeOfTask::kTBB;
  }
  explicit TimofeevNRadixBatcherTBB(const InType &in);

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

  void PrepAux(int &n, int &m, std::vector<int> &in, int &max_x, size_t &num_threads, size_t &n_thr);
  void HandleTBB(size_t &num_threads, size_t &n_n, size_t &m_m, int &max_x, std::vector<int> &reff,
                 std::vector<int> &r_in);
};

}  // namespace timofeev_n_radix_batcher_sort_threads
