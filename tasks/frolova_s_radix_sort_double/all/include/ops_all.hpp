#pragma once

#include <memory>
#include <vector>

#include "frolova_s_radix_sort_double/common/include/common.hpp"
#include "task/include/task.hpp"

namespace frolova_s_radix_sort_double {

class FrolovaSRadixSortDoubleALL : public BaseTask {
 public:
  using InType = std::vector<double>;
  using OutType = std::vector<double>;

  explicit FrolovaSRadixSortDoubleALL(const InType &in);

  bool ValidationImpl() override;
  bool PreProcessingImpl() override;
  bool RunImpl() override;
  bool PostProcessingImpl() override;
};

}  // namespace frolova_s_radix_sort_double
