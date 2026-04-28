#pragma once

#include "task/include/task.hpp"
#include "vinyaikina_e_multidimensional_integrals_simpson_method/common/include/common.hpp"

namespace vinyaikina_e_multidimensional_integrals_simpson_method {

class VinyaikinaEMultidimIntegrSimpsonTBB : public BaseTask {
 public:
  static constexpr ppc::task::TypeOfTask GetStaticTypeOfTask() {
    return ppc::task::TypeOfTask::kTBB;
  }
  explicit VinyaikinaEMultidimIntegrSimpsonTBB(const InType &in);

 private:
  bool ValidationImpl() override;
  bool PreProcessingImpl() override;
  bool RunImpl() override;
  bool PostProcessingImpl() override;

  double I_res_{0.0};
};

}  // namespace vinyaikina_e_multidimensional_integrals_simpson_method
