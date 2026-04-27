#pragma once

#include "savva_d_monte_carlo/common/include/common.hpp"
#include "task/include/task.hpp"

namespace savva_d_monte_carlo {

class SavvaDMonteCarloTBB : public BaseTask {
 public:
  static constexpr ppc::task::TypeOfTask GetStaticTypeOfTask() {
    return ppc::task::TypeOfTask::kTBB;
  }
  explicit SavvaDMonteCarloTBB(const InType &in);

 private:
  bool ValidationImpl() override;
  bool PreProcessingImpl() override;
  bool RunImpl() override;
  bool PostProcessingImpl() override;
};

}  // namespace savva_d_monte_carlo
