#pragma once

#include <vector>

#include "lifanov_k_simple_hoar_seq/common/include/common.hpp"
#include "task/include/task.hpp"

namespace lifanov_k_simple_hoar_seq {

class LifanovKSimpleHoarSEQ : public BaseTask {
public:
  static constexpr ppc::task::TypeOfTask GetStaticTypeOfTask() {
    return ppc::task::TypeOfTask::kSEQ;
  }
  explicit LifanovKSimpleHoarSEQ(const InType &in);

private:
  bool ValidationImpl() override;
  bool PreProcessingImpl() override;
  bool RunImpl() override;
  bool PostProcessingImpl() override;

  static void QuickSortHoare(std::vector<int> &arr, int low, int high);
  static int Partition(std::vector<int> &arr, int low, int high);
  static std::vector<int> Merge(const std::vector<int> &left,
                                const std::vector<int> &right);

  std::vector<int> data_;
};

} // namespace lifanov_k_simple_hoar_seq