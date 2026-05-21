#pragma once

#include <cstddef>
#include <vector>

#include "Nazarova_K_rad_sort_batcher_metod/common/include/common.hpp"
#include "task/include/task.hpp"

namespace nazarova_k_calc_integ_rectangles {

class NazarovaKCalcIntegRectanglesSTL : public BaseTask {
 public:
  static constexpr ppc::task::TypeOfTask GetStaticTypeOfTask() {
    return ppc::task::TypeOfTask::kSTL;
  }
  explicit NazarovaKCalcIntegRectanglesSTL(const InType &in);

 private:
  bool ValidationImpl() override;
  bool PreProcessingImpl() override;
  bool RunImpl() override;
  bool PostProcessingImpl() override;

  std::size_t dimension_{0U};
  std::vector<double> step_sizes_;
  double cell_volume_{0.0};
  std::size_t total_cells_{0U};
  double result_{0.0};
};

}  // namespace nazarova_k_calc_integ_rectangles
