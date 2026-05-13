#pragma once

#include <cstddef>
#include <vector>

#include "../../common/include/common.hpp"
#include "task/include/task.hpp"

namespace kutergin_v_multidimensional_integration_rect_method {

class RectMethodOMP : public BaseTask {
 public:
  static constexpr ppc::task::TypeOfTask GetStaticTypeOfTask() {
    return ppc::task::TypeOfTask::kOMP;
  }

  explicit RectMethodOMP(const InType &in);

 protected:
  bool ValidationImpl() override;
  bool PreProcessingImpl() override;
  bool RunImpl() override;
  bool PostProcessingImpl() override;

 private:
  double CalculateChunkSum(size_t start_idx, size_t count, const std::vector<double> &h);

  InType local_input_;
  double res_{0.0};
};

}  // namespace kutergin_v_multidimensional_integration_rect_method
