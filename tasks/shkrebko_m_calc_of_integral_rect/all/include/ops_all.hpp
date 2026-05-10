#pragma once

#include <utility>
#include <vector>

#include "shkrebko_m_calc_of_integral_rect/common/include/common.hpp"
#include "task/include/task.hpp"

namespace shkrebko_m_calc_of_integral_rect {

class ShkrebkoMCalcOfIntegralRectALL : public BaseTask {
 public:
  static constexpr ppc::task::TypeOfTask GetStaticTypeOfTask() {
    return ppc::task::TypeOfTask::kALL;
  }

  explicit ShkrebkoMCalcOfIntegralRectALL(const InType &in);

 protected:
  bool ValidationImpl() override;
  bool PreProcessingImpl() override;
  bool RunImpl() override;
  bool PostProcessingImpl() override;

 private:
  void BroadcastCommonData(int rank);
  void AssignMpiSlice(int rank, int size, double &local_left, double &local_right, int &local_steps, int &local_offset);

  [[nodiscard]] double ComputeSliceSum(double left, double right, int steps, const std::vector<double> &h_other,
                                       const std::vector<std::pair<double, double>> &limits_other,
                                       const std::vector<int> &n_steps_other) const;

  InType local_input_;
  double res_ = 0.0;
};

}  // namespace shkrebko_m_calc_of_integral_rect
