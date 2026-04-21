// ops_tbb.hpp
#pragma once

#include <tbb/tbb.h>

#include <vector>

#include "task/include/task.hpp"
#include "tsibareva_e_integral_calculate_trapezoid_method/common/include/common.hpp"

namespace tsibareva_e_integral_calculate_trapezoid_method {

class TsibarevaEIntegralCalculateTrapezoidMethodTBB : public BaseTask {
 public:
  static constexpr ppc::task::TypeOfTask GetStaticTypeOfTask() {
    return ppc::task::TypeOfTask::kTBB;
  }
  explicit TsibarevaEIntegralCalculateTrapezoidMethodTBB(const InType &in);

 private:
  bool ValidationImpl() override;
  bool PreProcessingImpl() override;
  bool RunImpl() override;
  bool PostProcessingImpl() override;

  static double ComputeSumNode(const tbb::blocked_range<int> &range, double init, const Integral &input,
                               const std::vector<double> &h, const std::vector<int> &sizes);
};

}  // namespace tsibareva_e_integral_calculate_trapezoid_method
