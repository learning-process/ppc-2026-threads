#pragma once

#include <cstddef>
#include <vector>

#include "Nazarova_K_calc_integ_rectangles/common/include/common.hpp"
#include "task/include/task.hpp"

namespace nazarova_k_calc_integ_rectangles {

class NazarovaKCalcIntegRectanglesOMP : public BaseTask {
 public:
  static constexpr ppc::task::TypeOfTask GetStaticTypeOfTask() {
    return ppc::task::TypeOfTask::kOMP;
  }
  explicit NazarovaKCalcIntegRectanglesOMP(const InType &in);

 private:
  bool ValidationImpl() override;
  bool PreProcessingImpl() override;
  bool RunImpl() override;
  bool PostProcessingImpl() override;

  [[nodiscard]] bool HasValidInput();

  std::size_t dimension_{0U};
  std::vector<double> step_sizes_;
  double cell_volume_{0.0};
  std::size_t total_cells_{0U};
  double result_{0.0};
};

}  // namespace nazarova_k_calc_integ_rectangles
