#pragma once

#include "task/include/task.hpp"
#include "tsibareva_e_integral_calculate_trapezoid_method/common/include/common.hpp"

namespace tsibareva_e_integral_calculate_trapezoid_method {

class TsibarevaEIntegralCalculateTrapezoidMethodSEQ : public ppc::task::Task<IntegralInput, double> {
 public:
 static constexpr ppc::task::TypeOfTask GetStaticTypeOfTask() {
    return ppc::task::TypeOfTask::kSEQ;
  }
  TsibarevaEIntegralCalculateTrapezoidMethodSEQ(const InType& in);  
  
  private:
  bool ValidationImpl() override;
  bool PreProcessingImpl() override;
  bool RunImpl() override;
  bool PostProcessingImpl() override;
};

}  // namespace tsibareva_e_integral_calculate_trapezoid_method