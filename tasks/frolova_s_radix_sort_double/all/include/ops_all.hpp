#pragma once

#include <vector>

#include "frolova_s_radix_sort_double/common/include/common.hpp"
#include "task/include/task.hpp"

namespace frolova_s_radix_sort_double {

class FrolovaSRadixSortDoubleALL : public ppc::core::Task {
 public:
  static constexpr ppc::task::TypeOfTask GetStaticTypeOfTask() {
    return ppc::task::TypeOfTask::kALL;
  }
  explicit FrolovaSRadixSortDoubleALL(const InType &in);

 private:
  bool ValidationImpl() override;
  bool PreProcessingImpl() override;
  bool RunImpl() override;
  bool PostProcessingImpl() override;

  std::vector<double> MergeSorted(const std::vector<double> &a, const std::vector<double> &b);
};

}  // namespace frolova_s_radix_sort_double
