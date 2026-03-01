#pragma once

#include "kutergin_a_multidim_trapezoid/common/include/common.hpp"
#include "task/include/task.hpp"

namespace kutergin_a_multidim_trapezoid {

class KuterginAMultidimTrapezoidSEQ : public BaseTask {
 public:
  static constexpr ppc::task::TypeOfTask GetStaticTypeOfTask() {
    return ppc::task::TypeOfTask::kSEQ;
  }

  explicit KuterginAMultidimTrapezoidSEQ(const InType &in);

 private:
  bool ValidationImpl() override;
  bool PreProcessingImpl() override;
  bool RunImpl() override;
  bool PostProcessingImpl() override;
};

}  // namespace kutergin_a_multidim_trapezoid
