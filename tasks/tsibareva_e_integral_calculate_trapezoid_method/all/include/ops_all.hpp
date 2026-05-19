#pragma once

#include <functional>
#include <vector>

#include "task/include/task.hpp"
#include "tsibareva_e_integral_calculate_trapezoid_method/common/include/common.hpp"

namespace tsibareva_e_integral_calculate_trapezoid_method {

class TsibarevaEIntegralCalculateTrapezoidMethodALL : public BaseTask {
 public:
  static constexpr ppc::task::TypeOfTask GetStaticTypeOfTask() {
    return ppc::task::TypeOfTask::kALL;
  }
  explicit TsibarevaEIntegralCalculateTrapezoidMethodALL(const InType &in);

 private:
  bool ValidationImpl() override;
  bool PreProcessingImpl() override;
  bool RunImpl() override;
  bool PostProcessingImpl() override;
  static double ComputePartialSum(int begin, int finish, const std::vector<double> &lo, const std::vector<double> &h,
                                  const std::vector<int> &sizes, const std::vector<int> &steps, int dim,
                                  const std::function<double(const std::vector<double> &)> &f);
};

}  // namespace tsibareva_e_integral_calculate_trapezoid_method
