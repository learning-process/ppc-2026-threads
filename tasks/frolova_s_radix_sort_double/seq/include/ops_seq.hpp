#pragma once

#include <vector>

#include "frolova_s_radix_sort_double/common/include/common.hpp"
#include "task/include/task.hpp"

namespace frolova_s_radix_sort_double {

class FrolovaSRadixSortDoubleSEQ : public BaseTask {
 public:
  static constexpr ppc::task::TypeOfTask GetStaticTypeOfTask() {
    return ppc::task::TypeOfTask::kSEQ;
  }
  explicit FrolovaSRadixSortDoubleSEQ(const InType &in);

 private:
  bool ValidationImpl() override;
  bool PreProcessingImpl() override;
  bool RunImpl() override;
  bool PostProcessingImpl() override;
};

}  // namespace frolova_s_radix_sort_double