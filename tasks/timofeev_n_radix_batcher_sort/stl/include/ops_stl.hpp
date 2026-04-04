#pragma once

#include <cstddef>
#include <vector>

#include "task/include/task.hpp"
#include "timofeev_n_radix_batcher_sort/common/include/common.hpp"

namespace timofeev_n_radix_batcher_sort_threads {

class TimofeevNRadixBatcherSTL : public BaseTask {
 public:
  static constexpr ppc::task::TypeOfTask GetStaticTypeOfTask() {
    return ppc::task::TypeOfTask::kSTL;
  }
  explicit TimofeevNRadixBatcherSTL(const InType &in);

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
  static void BubbleSortAux(size_t &num_threads, int &max_x, std::vector<int> &r_in, size_t &piece,
                            std::vector<std::thread> &threads);
  static void OddMergeAux(std::vector<int> &r_in, size_t &piece, std::vector<std::thread> &threads, size_t &n_n);
  static void HandleSTL(size_t &num_threads, size_t &n_n, size_t &m_m, int &max_x, std::vector<int> &reff,
                        std::vector<int> &r_in);
};

}  // namespace timofeev_n_radix_batcher_sort_threads
