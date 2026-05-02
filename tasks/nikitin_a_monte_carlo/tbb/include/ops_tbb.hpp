#pragma once

#include "nikitin_a_monte_carlo/common/include/common.hpp"
#include "task/include/task.hpp"

namespace nikitin_a_monte_carlo {

class NikitinAMonteCarloTBB : public BaseTask {
 public:
  static constexpr ppc::task::TypeOfTask GetStaticTypeOfTask() {
    return ppc::task::TypeOfTask::kTBB;
  }
  explicit NikitinAMonteCarloTBB(const InType &in);

 private:
  bool ValidationImpl() override;
  bool PreProcessingImpl() override;
  bool RunImpl() override;
  bool PostProcessingImpl() override;
};

}  // namespace nikitin_a_monte_carlo
