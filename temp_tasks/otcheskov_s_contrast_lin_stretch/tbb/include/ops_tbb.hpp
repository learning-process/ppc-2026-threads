#pragma once

#include <tbb/tbb.h>

#include <cstdint>

#include "otcheskov_s_contrast_lin_stretch/common/include/common.hpp"
#include "task/include/task.hpp"

namespace otcheskov_s_contrast_lin_stretch {

class OtcheskovSContrastLinStretchTBB : public BaseTask {
 public:
  static constexpr ppc::task::TypeOfTask GetStaticTypeOfTask() {
    return ppc::task::TypeOfTask::kTBB;
  }
  explicit OtcheskovSContrastLinStretchTBB(const InType &in);

 private:
  bool ValidationImpl() override;
  bool PreProcessingImpl() override;
  bool RunImpl() override;
  bool PostProcessingImpl() override;

  struct MinMax {
    uint8_t min{255};
    uint8_t max{0};
  };

  static MinMax ComputeMinMax(const InType &input, tbb::task_arena &arena);
  static void CopyInput(const InType &input, OutType &output, tbb::task_arena &arena);
  static void LinearStretch(const InType &input, OutType &output, int min_i, int range, tbb::task_arena &arena);
};

}  // namespace otcheskov_s_contrast_lin_stretch
