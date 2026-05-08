#pragma once

#include <cstddef>
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
  void BroadcastInputData();
  [[nodiscard]] std::size_t SelectSplitDimension() const;
  [[nodiscard]] double ComputeSliceSum(std::size_t fixed_dim, std::size_t fixed_idx,
                                       const std::vector<double> &h) const;

  InType local_input_;
  double res_ = 0.0;
  int rank_ = 0;
  int world_size_ = 1;
  bool use_mpi_ = false;
};

}  // namespace shkrebko_m_calc_of_integral_rect
