#pragma once

#include <utility>
#include <vector>

#include "task/include/task.hpp"
#include "tochilin_e_hoar_sort_sim_mer_seq/common/include/common.hpp"

namespace tochilin_e_hoar_sort_sim_mer_seq {

class TochilinEHoarSortSimMerSEQ : public BaseTask {
public:
  static constexpr ppc::task::TypeOfTask GetStaticTypeOfTask() {
    return ppc::task::TypeOfTask::kSEQ;
  }
  explicit TochilinEHoarSortSimMerSEQ(const InType &in);

private:
  bool ValidationImpl() override;
  bool PreProcessingImpl() override;
  bool RunImpl() override;
  bool PostProcessingImpl() override;

  static void QuickSort(std::vector<int> &arr, int low, int high);
  static std::pair<int, int> Partition(std::vector<int> &arr, int l, int r);
  static std::vector<int> MergeSortedVectors(const std::vector<int> &a,
                                             const std::vector<int> &b);
};

} // namespace tochilin_e_hoar_sort_sim_mer_seq
