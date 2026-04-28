#pragma once

#include <cstdint>

#include "otcheskov_s_contrast_lin_stretch/common/include/common.hpp"
#include "task/include/task.hpp"

namespace otcheskov_s_contrast_lin_stretch {

class OtcheskovSContrastLinStretchSTL : public BaseTask {
 public:
  static constexpr ppc::task::TypeOfTask GetStaticTypeOfTask() {
    return ppc::task::TypeOfTask::kSTL;
  }
  explicit OtcheskovSContrastLinStretchSTL(const InType &in);

 private:
  bool ValidationImpl() override;
  bool PreProcessingImpl() override;
  bool RunImpl() override;
  bool PostProcessingImpl() override;

  struct MinMax {
    uint8_t min{255};
    uint8_t max{0};
  };

  static MinMax ComputeMinMax(const InType &input);
  static void CopyInput(const InType &input, OutType &output);
  static void LinearStretch(const InType &input, OutType &output, int min_i, int range);
};

}  // namespace otcheskov_s_contrast_lin_stretch
