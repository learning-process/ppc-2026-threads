#pragma once

#include <vector>

#include "chernov_t_radix_sort/common/include/common.hpp"
#include "task/include/task.hpp"

namespace chernov_t_radix_sort {

class ChernovTRadixSortOMP : public BaseTask {
 public:
  static constexpr ppc::task::TypeOfTask GetStaticTypeOfTask() {
    return ppc::task::TypeOfTask::kOMP;
  }
  explicit ChernovTRadixSortOMP(const InType &in);

 private:
  bool ValidationImpl() override;
  bool PreProcessingImpl() override;
  bool RunImpl() override;
  bool PostProcessingImpl() override;

  static void RadixSortLSD(std::vector<int> &data);
};

}  // namespace chernov_t_radix_sort
