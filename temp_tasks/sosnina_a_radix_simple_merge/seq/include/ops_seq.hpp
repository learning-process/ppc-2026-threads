#pragma once

#include "sosnina_a_radix_simple_merge/common/include/common.hpp"
#include "task/include/task.hpp"

namespace sosnina_a_radix_simple_merge {

class SosninaATestTaskSEQ : public BaseTask {
 public:
  static constexpr ppc::task::TypeOfTask GetStaticTypeOfTask() {
    return ppc::task::TypeOfTask::kSEQ;
  }
  explicit SosninaATestTaskSEQ(const InType &in);

 private:
  bool ValidationImpl() override;
  bool PreProcessingImpl() override;
  bool RunImpl() override;
  bool PostProcessingImpl() override;
};

}  // namespace sosnina_a_radix_simple_merge
