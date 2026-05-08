#pragma once
#include <cstddef>
#include <vector>

#include "chernykh_s_trapezoidal_integration/common/include/common.hpp"
#include "task/include/task.hpp"
namespace chernykh_s_trapezoidal_integration {

class ChernykhSTrapezoidalIntegrationTBB : public BaseTask {
 public:
  static constexpr ppc::task::TypeOfTask GetStaticTypeOfTask() {
    return ppc::task::TypeOfTask::kTBB;
  }
  explicit ChernykhSTrapezoidalIntegrationTBB(const InType &in);

 private:
  static double CalculatePointAndWeight(const IntegrationInType &input, const std::vector<std::size_t> &counters,
                                        std::vector<double> &point);

  bool ValidationImpl() override;
  bool PreProcessingImpl() override;
  bool RunImpl() override;
  bool PostProcessingImpl() override;
};

}  // namespace chernykh_s_trapezoidal_integration
