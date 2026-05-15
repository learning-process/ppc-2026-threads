#pragma once

#include <vector>

#include "task/include/task.hpp"

namespace frolova_s_radix_sort_double {

class FrolovaSRadixSortDoubleALL : public ppc::task::BaseTask {
 public:
  using InType = std::vector<double>;
  using OutType = std::vector<double>;

  static constexpr ppc::task::TypeOfTask GetStaticTypeOfTask() {
    return ppc::task::TypeOfTask::kALL;
  }

  explicit FrolovaSRadixSortDoubleALL(const InType &in);

  bool ValidationImpl() override;
  bool PreProcessingImpl() override;
  bool RunImpl() override;
  bool PostProcessingImpl() override;

 private:
  InType input_;
  OutType output_;
};

}  // namespace frolova_s_radix_sort_double
