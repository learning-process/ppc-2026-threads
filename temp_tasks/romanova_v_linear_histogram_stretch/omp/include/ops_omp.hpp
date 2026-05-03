#pragma once

#include "romanova_v_linear_histogram_stretch/common/include/common.hpp"
#include "task/include/task.hpp"

namespace romanova_v_linear_histogram_stretch_threads {

class RomanovaVLinHistogramStretchOMP : public BaseTask {
 public:
  static constexpr ppc::task::TypeOfTask GetStaticTypeOfTask() {
    return ppc::task::TypeOfTask::kOMP;
  }
  explicit RomanovaVLinHistogramStretchOMP(const InType &in);

 private:
  bool ValidationImpl() override;
  bool PreProcessingImpl() override;
  bool RunImpl() override;
  bool PostProcessingImpl() override;
};

}  // namespace romanova_v_linear_histogram_stretch_threads
