#pragma once

#include <utility>
#include <vector>

#include "task/include/task.hpp"
#include "tochilin_e_hoar_sort_sim_mer/common/include/common.hpp"

namespace tochilin_e_hoar_sort_sim_mer {

class TochilinEHoarSortSimMerOMP : public BaseTask {
 public:
  static constexpr ppc::task::TypeOfTask GetStaticTypeOfTask() {
    return ppc::task::TypeOfTask::kOMP;
  }
  explicit TochilinEHoarSortSimMerOMP(const InType &in);

 private:
  bool ValidationImpl() override;
  bool PreProcessingImpl() override;
  bool RunImpl() override;
  bool PostProcessingImpl() override;

  static void QuickSortOMP(std::vector<int> &arr, int low, int high, int depth_limit);
  static std::pair<int, int> Partition(std::vector<int> &arr, int l, int r);
  static std::vector<int> MergeSortedVectors(const std::vector<int> &a, const std::vector<int> &b);
};

}  // namespace tochilin_e_hoar_sort_sim_mer
