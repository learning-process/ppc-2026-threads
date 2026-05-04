#pragma once

#include "samoylenko_i_integral_trapezoid/common/include/common.hpp"
#include "task/include/task.hpp"

namespace samoylenko_i_integral_trapezoid {

class SamoylenkoIIntegralTrapezoidOMP : public BaseTask {
 public:
  static constexpr ppc::task::TypeOfTask GetStaticTypeOfTask() {
    return ppc::task::TypeOfTask::kOMP;
  }

  explicit SamoylenkoIIntegralTrapezoidOMP(const InType &in);

 private:
  bool ValidationImpl() override;
  bool PreProcessingImpl() override;
  bool RunImpl() override;
  bool PostProcessingImpl() override;
};

}  // namespace samoylenko_i_integral_trapezoid
