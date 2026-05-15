#pragma once

#include <vector>

#include "chaschin_vladimir_linear_image_filtration_seq/common/include/common.hpp"
#include "task/include/task.hpp"

namespace chaschin_v_linear_image_filtration_all {

class ChaschinVLinearFiltrationALL : public chaschin_v_linear_image_filtration_seq::BaseTask {
 public:
  static constexpr ppc::task::TypeOfTask GetStaticTypeOfTask() {
    return ppc::task::TypeOfTask::kALL;
  }
  explicit ChaschinVLinearFiltrationALL(const chaschin_v_linear_image_filtration_seq::InType &in);

 private:
  bool ValidationImpl() override;
  bool PreProcessingImpl() override;
  bool RunImpl() override;
  bool PostProcessingImpl() override;
};

inline float HorizontalFilterAtALL(const std::vector<float> &img, int n, int x, int y);
inline float VerticalFilterAtALL(const std::vector<float> &temp, int n, int m, int x, int y);

}  // namespace chaschin_v_linear_image_filtration_all
