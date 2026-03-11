#pragma once

#include "pikhotskiy_r_vertical_gauss_filter/common/include/common.hpp"
#include "task/include/task.hpp"

namespace pikhotskiy_r_vertical_gauss_filter {

class PikhotskiyRVerticalGaussFilterSEQ : public BaseTask {
 public:
  static constexpr ppc::task::TypeOfTask GetStaticTypeOfTask() {
    return ppc::task::TypeOfTask::kSEQ;
  }
  explicit PikhotskiyRVerticalGaussFilterSEQ(const InType &in);

 private:
  bool ValidationImpl() override;
  bool PreProcessingImpl() override;
  bool RunImpl() override;
  bool PostProcessingImpl() override;
};

}  // namespace pikhotskiy_r_vertical_gauss_filter
// Пустая строка в конце
