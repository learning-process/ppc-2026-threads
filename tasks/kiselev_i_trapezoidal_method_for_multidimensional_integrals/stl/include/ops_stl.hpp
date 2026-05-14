#pragma once

#include <vector>

#include "kiselev_i_trapezoidal_method_for_multidimensional_integrals/common/include/common.hpp"
#include "task/include/task.hpp"

namespace kiselev_i_trapezoidal_method_for_multidimensional_integrals {

class KiselevITestTaskSTL : public BaseTask {
 public:
  static constexpr ppc::task::TypeOfTask GetStaticTypeOfTask() {
    return ppc::task::TypeOfTask::kSTL;
  }
  explicit KiselevITestTaskSTL(const InType &in);

  static double FunctionTypeChoose(int type_x, double x, double y);

 private:
  bool ValidationImpl() override;
  bool PreProcessingImpl() override;
  bool RunImpl() override;
  bool PostProcessingImpl() override;
  double ComputeIntegral(const std::vector<int> &steps);
};

}  // namespace kiselev_i_trapezoidal_method_for_multidimensional_integrals
