#pragma once

#include "kruglova_a_conjugate_gradient_sle/common/include/common.hpp"
#include "task/include/task.hpp"

namespace kruglova_a_conjugate_gradient_sle {

class KruglovaAConjGradSleSTL : public BaseTask {
 public:
  static constexpr ppc::task::TypeOfTask GetStaticTypeOfTask() {
    return ppc::task::TypeOfTask::kSTL;
  }
  explicit KruglovaAConjGradSleSTL(const InType &in);

 private:
  bool ValidationImpl() override;
  bool PreProcessingImpl() override;
  bool RunImpl() override;
  bool PostProcessingImpl() override;
};

}  // namespace kruglova_a_conjugate_gradient_sle
