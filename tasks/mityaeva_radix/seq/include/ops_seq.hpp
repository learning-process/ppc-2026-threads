#pragma once

#include "mityaeva_radix/common/include/common.hpp"
#include "task/include/task.hpp"

namespace mityaeva_radix {

class MityaevaRadixSeq : public BaseTask {
 public:
  static constexpr ppc::task::TypeOfTask GetStaticTypeOfTask() {
    return ppc::task::TypeOfTask::kSEQ;
  }
  explicit MityaevaRadixSeq(const InType &in);

 private:
  bool ValidationImpl() override;
  bool PreProcessingImpl() override;
  bool RunImpl() override;
  bool PostProcessingImpl() override;
};

}  // namespace mityaeva_radix
