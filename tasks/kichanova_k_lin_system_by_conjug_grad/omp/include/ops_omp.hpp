#pragma once

#include "kichanova_k_lin_system_by_conjug_grad/common/include/common.hpp"
#include "task/include/task.hpp"

namespace kichanova_k_lin_system_by_conjug_grad {

class KichanovaKLinSystemByConjugGradOMP : public BaseTask {
 public:
  static constexpr ppc::task::TypeOfTask GetStaticTypeOfTask() {
    return ppc::task::TypeOfTask::kOMP;
  }
  explicit KichanovaKLinSystemByConjugGradOMP(const InType &in);

 private:
  bool ValidationImpl() override;
  bool PreProcessingImpl() override;
  bool RunImpl() override;
  bool PostProcessingImpl() override;
};

}  // namespace kichanova_k_lin_system_by_conjug_grad
