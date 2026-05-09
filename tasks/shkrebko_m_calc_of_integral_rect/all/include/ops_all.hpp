#pragma once

#include <cstddef>
#include <functional>
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
  void BroadcastInputData(int rank, std::size_t &dims);
  std::size_t SelectSplitDimension() const;
  bool ComputeSliceSum(std::size_t fixed_dim, std::size_t fixed_idx, const std::vector<double> &h,
                       double &slice_sum) const;

  InType local_input_;
  double res_ = 0.0;
};

}  // namespace shkrebko_m_calc_of_integral_rect
