#pragma once

#include <cstddef>
#include <cstdint>
#include <vector>

#include "chernykh_s_trapezoidal_integration/common/include/common.hpp"
#include "task/include/task.hpp"

namespace chernykh_s_trapezoidal_integration {

class ChernykhSTrapezoidalIntegrationALL : public BaseTask {
 public:
  static constexpr ppc::task::TypeOfTask GetStaticTypeOfTask() {
    return ppc::task::TypeOfTask::kALL;
  }
  explicit ChernykhSTrapezoidalIntegrationALL(const InType &in);

 private:
  static double CalculatePointAndWeight(const IntegrationInType &input, const std::vector<std::size_t> &counters,
                                        std::vector<double> &point);
  static double OnProcessCalculate(const IntegrationInType &input, std::size_t dims, int64_t start, int64_t end);
  static int DetectFunctionId(const IntegrationInType &input, std::size_t dims);

  bool ValidationImpl() override;
  bool PreProcessingImpl() override;
  bool RunImpl() override;
  bool PostProcessingImpl() override;
};

}  // namespace chernykh_s_trapezoidal_integration
