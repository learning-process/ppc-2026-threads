#pragma once

#include <vector>

#include "task/include/task.hpp"
#include "tsibareva_e_integral_calculate_trapezoid_method/common/include/common.hpp"

namespace tsibareva_e_integral_calculate_trapezoid_method {

class TsibarevaEIntegralCalculateTrapezoidMethodSEQ : public ppc::task::Task<Integral, double> {
 public:
  static constexpr ppc::task::TypeOfTask GetStaticTypeOfTask() {
    return ppc::task::TypeOfTask::kSEQ;
  }
  explicit TsibarevaEIntegralCalculateTrapezoidMethodSEQ(const InType &in);

 private:
  bool ValidationImpl() override;
  bool PreProcessingImpl() override;
  bool RunImpl() override;
  bool PostProcessingImpl() override;

  std::vector<double> ComputePoint(const std::vector<int> &indexes, const std::vector<double> &h, int dim);
  bool IterateGridPoints(std::vector<int> &indexes, int dim);
};

}  // namespace tsibareva_e_integral_calculate_trapezoid_method
