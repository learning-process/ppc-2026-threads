#pragma once

#include <cstddef>
#include <thread>
#include <vector>

#include "task/include/task.hpp"
#include "timofeev_n_radix_batcher_sort/common/include/common.hpp"

namespace timofeev_n_radix_batcher_sort_threads {

class TimofeevNRadixBatcherALL : public BaseTask {
 public:
  static constexpr ppc::task::TypeOfTask GetStaticTypeOfTask() {
    return ppc::task::TypeOfTask::kALL;
  }
  explicit TimofeevNRadixBatcherALL(const InType &in);

 private:
  bool ValidationImpl() override;
  bool PreProcessingImpl() override;
  bool RunImpl() override;
  bool PostProcessingImpl() override;

  static void CompExch(int &a, int &b, int digit);
  static void BubbleSort(std::vector<int> &arr, int digit, int left, int right);
  static void ComparR(int &a, int &b);
  static void OddEvenMerge(std::vector<int> &arr, int lft, int n);
  static void BubbleSortAux(size_t &num_threads, int &max_x, std::vector<int> &r_in, size_t &piece,
                            std::vector<std::thread> &threads);
  static void OddMergeAux(std::vector<int> &r_in, size_t &piece, std::vector<std::thread> &threads, size_t &n_n);
  static void ProcessLocalArray(std::vector<int> &local_arr, size_t num_threads);
  static void ProcessLocalArrayWOSort(std::vector<int> &local_arr, size_t num_threads, size_t &elements_per_process);
  static void PrepAux(int &n, int &m, std::vector<int> &in, int &max_x, size_t &num_threads, size_t &n_thr,
                      size_t &number_procs);
  void HandleZero(std::vector<int> &global_array, size_t &total_elements, size_t &total_elements_primal,
                  int &num_processes, std::vector<int> &local_array, int &maxxx, size_t &num_threads_per_process,
                  size_t &elements_per_process);
};

}  // namespace timofeev_n_radix_batcher_sort_threads
