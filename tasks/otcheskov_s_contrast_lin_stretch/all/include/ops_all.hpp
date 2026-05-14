#pragma once

#include <cstdint>

#include "otcheskov_s_contrast_lin_stretch/common/include/common.hpp"
#include "task/include/task.hpp"

namespace otcheskov_s_contrast_lin_stretch {

class OtcheskovSContrastLinStretchALL : public BaseTask {
 public:
  static constexpr ppc::task::TypeOfTask GetStaticTypeOfTask() {
    return ppc::task::TypeOfTask::kALL;
  }
  explicit OtcheskovSContrastLinStretchALL(const InType &in);

 private:
  int rank_{};
  int size_{};
  bool is_valid_{};
  struct MinMax {
    uint8_t min{255};
    uint8_t max{0};
  };

  bool ValidationImpl() override;
  bool PreProcessingImpl() override;
  bool RunImpl() override;
  bool PostProcessingImpl() override;

  static MinMax ComputeMinMax(const InType &input);
  static void CopyInput(const InType &input, OutType &output);
  static void LinearStretch(const InType &input, OutType &output, int min_i, int range);
};

}  // namespace otcheskov_s_contrast_lin_stretch
