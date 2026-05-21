#pragma once

#include <cstddef>
#include <vector>

#include "chernov_t_radix_sort/common/include/common.hpp"
#include "task/include/task.hpp"

namespace chernov_t_radix_sort {

class ChernovTRadixSortALL : public BaseTask {
 public:
  static constexpr ppc::task::TypeOfTask GetStaticTypeOfTask() {
    return ppc::task::TypeOfTask::kALL;
  }
  explicit ChernovTRadixSortALL(const InType &in);

 private:
  bool ValidationImpl() override;
  bool PreProcessingImpl() override;
  bool RunImpl() override;
  bool PostProcessingImpl() override;

  static void RadixSortLSDParallelOMP(std::vector<int> &data);
  static void SimpleMerge(const std::vector<int> &left, const std::vector<int> &right, std::vector<int> &result);

  static void ComputeChunkSizes(int num_processes, size_t total_elements, std::vector<int> &recv_counts,
                                std::vector<int> &displs);
  static void MergeChunksOnRank0(const std::vector<int> &global_result, const std::vector<int> &recv_counts,
                                 const std::vector<int> &displs, std::vector<int> &output);
};

}  // namespace chernov_t_radix_sort
