#pragma once

#include <vector>

#include "sabirov_s_monte_carlo_seq/common/include/common.hpp"
#include "task/include/task.hpp"

namespace sabirov_s_monte_carlo_seq {

class SabirovSMonteCarloSEQ : public BaseTask {
 public:
  static constexpr ppc::task::TypeOfTask GetStaticTypeOfTask() {
    return ppc::task::TypeOfTask::kSEQ;
  }
  explicit SabirovSMonteCarloSEQ(const InType &in);

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

}  // namespace sabirov_s_monte_carlo_seq
