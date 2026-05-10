#pragma once

#include <cstddef>
#include <functional>
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
  void DistributeData(int rank, size_t &dims);

  static double ComputeChunkSum(size_t start_idx, size_t end_idx, const std::vector<double> &h,
                                const std::vector<std::pair<double, double>> &limits, const std::vector<int> &n_steps,
                                const std::function<double(const std::vector<double> &)> &func);

  InType local_input_;
  double res_ = 0.0;
};

}  // namespace shkrebko_m_calc_of_integral_rect
