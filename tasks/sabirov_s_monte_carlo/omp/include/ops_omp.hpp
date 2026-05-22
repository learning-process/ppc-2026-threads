#pragma once

#include <vector>

#include "sabirov_s_monte_carlo/common/include/common.hpp"
#include "task/include/task.hpp"

namespace sabirov_s_monte_carlo {

class SabirovSMonteCarloOMP : public BaseTask {
 public:
  static constexpr ppc::task::TypeOfTask GetStaticTypeOfTask() {
    return ppc::task::TypeOfTask::kOMP;
  }
  explicit SabirovSMonteCarloOMP(const InType &in);

 private:
  bool ValidationImpl() override;
  bool PreProcessingImpl() override;
  bool RunImpl() override;
  bool PostProcessingImpl() override;

  std::vector<double> lower_;
  std::vector<double> upper_;
  int num_samples_{};
  FuncType func_type_{};
};

}  // namespace sabirov_s_monte_carlo
