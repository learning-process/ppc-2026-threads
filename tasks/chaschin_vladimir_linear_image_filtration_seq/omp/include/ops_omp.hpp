#pragma once

#include <vector>

#include "chaschin_vladimir_linear_image_filtration_seq/common/include/common.hpp"
#include "task/include/task.hpp"

namespace chaschin_v_linear_image_filtration_omp {

class ChaschinVLinearFiltrationOMP : public chaschin_v_linear_image_filtration_seq::BaseTask {
 public:
  static constexpr ppc::task::TypeOfTask GetStaticTypeOfTask() {
    return ppc::task::TypeOfTask::kOMP;
  }
  explicit ChaschinVLinearFiltrationOMP(const chaschin_v_linear_image_filtration_seq::InType &in);

 private:
  bool ValidationImpl() override;
  bool PreProcessingImpl() override;
  bool RunImpl() override;
  bool PostProcessingImpl() override;
};

inline float HorizontalFilterAt(const std::vector<float> &img, int n, int x, int y);
inline float VerticalFilterAt(const std::vector<float> &temp, int n, int m, int x, int y);

}  // namespace chaschin_v_linear_image_filtration_omp
