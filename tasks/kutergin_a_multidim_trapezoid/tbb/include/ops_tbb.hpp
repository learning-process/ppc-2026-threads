#pragma once

#include "kutergin_a_multidim_trapezoid/common/include/common.hpp"
#include "task/include/task.hpp"

namespace kutergin_a_multidim_trapezoid {

class KuterginAMultidimTrapezoidTBB : public BaseTask {
 public:
  static constexpr ppc::task::TypeOfTask GetStaticTypeOfTask() {
    return ppc::task::TypeOfTask::kTBB;
  }

  explicit KuterginAMultidimTrapezoidTBB(const InType &in);

 private:
  bool ValidationImpl() override;
  bool PreProcessingImpl() override;
  bool RunImpl() override;
  bool PostProcessingImpl() override;

  InType local_input_;  // Храним копию входных данных
  double res_{0.0};     // Храним результат
};

}  // namespace kutergin_a_multidim_trapezoid
