#pragma once

#include <vector>

#include "chernov_t_radix_sort/common/include/common.hpp"
#include "task/include/task.hpp"

namespace chernov_t_radix_sort {

class ChernovTRadixSortTBB : public BaseTask {
 public:
  static constexpr ppc::task::TypeOfTask GetStaticTypeOfTask() {
    return ppc::task::TypeOfTask::kTBB;
  }
  explicit ChernovTRadixSortTBB(const InType &in);

 private:
  bool ValidationImpl() override;
  bool PreProcessingImpl() override;
  bool RunImpl() override;
  bool PostProcessingImpl() override;

  static void RadixSortLSD(std::vector<int> &data);
  static void SimpleMerge(const std::vector<int> &left, const std::vector<int> &right, std::vector<int> &result);
  static void ComputePrefixSums(std::vector<int> &count);
};

}  // namespace chernov_t_radix_sort
