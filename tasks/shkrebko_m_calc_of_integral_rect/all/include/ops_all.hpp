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
  void BroadcastInputData(int rank, std::size_t &dims);
  [[nodiscard]] std::size_t SelectSplitDimension() const;
  bool ComputeSliceSum(std::size_t fixed_dim, std::size_t fixed_idx, const std::vector<double> &h,
                       double &slice_sum) const;

  // Новые вспомогательные методы для снижения когнитивной сложности RunImpl
  [[nodiscard]] double ComputeCellVolume(const std::vector<double> &h) const;
  [[nodiscard]] std::vector<std::size_t> DistributeSlices(int rank, int size, std::size_t split_steps) const;
  [[nodiscard]] std::pair<double, bool> ComputeLocalSum(const std::vector<std::size_t> &local_slices,
                                                        const std::vector<double> &h, std::size_t split_dim) const;
  bool FinalizeResult(double local_sum, double cell_volume, bool local_ok, bool is_mpi, int rank);

  InType local_input_;
  double res_ = 0.0;
};

}  // namespace shkrebko_m_calc_of_integral_rect
