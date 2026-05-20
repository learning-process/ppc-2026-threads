#pragma once

#include <cstddef>
#include <utility>
#include <vector>

#include "task/include/task.hpp"
#include "tochilin_e_hoar_sort_sim_mer/common/include/common.hpp"

namespace tochilin_e_hoar_sort_sim_mer {

class TochilinEHoarSortSimMerSTL : public BaseTask {
 public:
  static constexpr ppc::task::TypeOfTask GetStaticTypeOfTask() {
    return ppc::task::TypeOfTask::kSTL;
  }
  explicit TochilinEHoarSortSimMerSTL(const InType &in);

 private:
  bool ValidationImpl() override;
  bool PreProcessingImpl() override;
  bool RunImpl() override;
  bool PostProcessingImpl() override;

  static void QuickSortSequential(std::vector<int> &arr, int low, int high);
  static void QuickSortParallel(std::vector<int> &arr, int low, int high, int serial_cutoff, int worker_count);
  static std::pair<int, int> Partition(std::vector<int> &arr, int l, int r);
  static bool ProcessRange(std::vector<int> &arr, std::pair<int, int> range, int serial_cutoff,
                           std::vector<std::pair<int, int>> &next_ranges);
  static int ResolveSerialCutoff(std::size_t size);
  static std::vector<int> MergeSortedVectors(const std::vector<int> &a, const std::vector<int> &b);
};

}  // namespace tochilin_e_hoar_sort_sim_mer
