#pragma once

#include <vector>

#include "shilin_n_monte_carlo_integration/common/include/common.hpp"
#include "task/include/task.hpp"

namespace shilin_n_monte_carlo_integration {

class ShilinNMonteCarloIntegrationSEQ : public BaseTask {
 public:
  static constexpr ppc::task::TypeOfTask GetStaticTypeOfTask() {
    return ppc::task::TypeOfTask::kSEQ;
  }
  explicit ShilinNMonteCarloIntegrationSEQ(const InType &in);

 private:
  bool ValidationImpl() override;
  bool PreProcessingImpl() override;
  bool RunImpl() override;
  bool PostProcessingImpl() override;

  std::vector<double> lower_bounds_;
  std::vector<double> upper_bounds_;
  int num_points_{};
  FuncType func_type_{};
};

}  // namespace shilin_n_monte_carlo_integration
