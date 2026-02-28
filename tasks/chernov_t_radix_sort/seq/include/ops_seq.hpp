#pragma once

#include "chernov_t_radix_sort/common/include/common.hpp"
#include "task/include/task.hpp"

namespace chernov_t_radix_sort {

class ChernovTRadixSortSEQ : public BaseTask {
 public:
  static constexpr ppc::task::TypeOfTask GetStaticTypeOfTask() {
    return ppc::task::TypeOfTask::kSEQ;
  }
  explicit ChernovTRadixSortSEQ(const InType &in);

 private:
  bool ValidationImpl() override;
  bool PreProcessingImpl() override;
  bool RunImpl() override;
  bool PostProcessingImpl() override;

  void RadixSortLSD(std::vector<int>& data);
  void SimpleMerge(const std::vector<int>& left,
                 const std::vector<int>& right,
                 std::vector<int>& result);
};

}  // namespace chernov_t_radix_sort
