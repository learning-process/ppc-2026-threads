#pragma once

#include "shkrebko_m_calc_of_integral_rect/common/include/common.hpp"
#include "task/include/task.hpp"

namespace shkrebko_m_calc_of_integral_rect {

class ShkrebkoMCalcOfIntegralRectSEQ : public BaseTask {
 public:
  static constexpr ppc::task::TypeOfTask GetStaticTypeOfTask() {
    return ppc::task::TypeOfTask::kSEQ;
  }

  explicit ShkrebkoMCalcOfIntegralRectSEQ(const InType &in);

 protected:
  bool ValidationImpl() override;
  bool PreProcessingImpl() override;
  bool RunImpl() override;
  bool PostProcessingImpl() override;

 private:
  InType local_input_;
  double res_;
};

}  // namespace shkrebko_m_calc_of_integral_rect
