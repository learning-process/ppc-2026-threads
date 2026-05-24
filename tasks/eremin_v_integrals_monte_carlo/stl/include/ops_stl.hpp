#pragma once

#include "eremin_v_integrals_monte_carlo/common/include/common.hpp"
#include "task/include/task.hpp"

namespace eremin_v_integrals_monte_carlo {

class EreminVIntegralsMonteCarloSTL : public BaseTask {
 public:
  static constexpr ppc::task::TypeOfTask GetStaticTypeOfTask() {
    return ppc::task::TypeOfTask::kSTL;
  }
  explicit EreminVIntegralsMonteCarloSTL(const InType &in);

 private:
  bool ValidationImpl() override;
  bool PreProcessingImpl() override;
  bool RunImpl() override;
  bool PostProcessingImpl() override;
};

}  // namespace eremin_v_integrals_monte_carlo
