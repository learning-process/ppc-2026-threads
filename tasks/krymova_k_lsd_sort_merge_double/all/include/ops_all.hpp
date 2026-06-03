#pragma once

#include <cstdint>
#include <vector>

#include "krymova_k_lsd_sort_merge_double/common/include/common.hpp"
#include "task/include/task.hpp"

namespace krymova_k_lsd_sort_merge_double {

class KrymovaKLsdSortMergeDoubleALL : public BaseTask {
 public:
  static constexpr ppc::task::TypeOfTask GetStaticTypeOfTask() {
    return ppc::task::TypeOfTask::kALL;
  }
  explicit KrymovaKLsdSortMergeDoubleALL(const InType &in);

 private:
  bool ValidationImpl() override;
  bool PreProcessingImpl() override;
  bool RunImpl() override;
  bool PostProcessingImpl() override;

  static void LSDSort(double *arr, int size);
  static std::vector<double> SimpleMerge(const std::vector<double> &a, const std::vector<double> &b);

  static uint64_t DoubleToULL(double d);
  static double ULLToDouble(uint64_t ull);

  bool RunSmallDataset(int total_size);
  static void ComputeDistribution(int total_size, int size_comm, std::vector<int> &send_counts,
                                  std::vector<int> &offsets);
  void ScatterData(int rank, const std::vector<int> &send_counts, const std::vector<int> &offsets,
                   std::vector<double> &local_data);
  void GatherResults(int rank, int size_comm, const std::vector<int> &send_counts, std::vector<double> &local_data);
  void BroadcastResult(int rank);
};

}  // namespace krymova_k_lsd_sort_merge_double
