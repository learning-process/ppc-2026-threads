#pragma once

#include "mityaeva_radix/common/include/common.hpp"
#include "task/include/task.hpp"

namespace mityaeva_radix {

class MityaevaRadixTbb : public BaseTask {
 public:
  static constexpr ppc::task::TypeOfTask GetStaticTypeOfTask() {
    return ppc::task::TypeOfTask::kTBB;
  }
  explicit MityaevaRadixTbb(const InType &in);

 private:
  bool ValidationImpl() override;
  bool PreProcessingImpl() override;
  bool RunImpl() override;
  bool PostProcessingImpl() override;
};

}  // namespace mityaeva_radix
