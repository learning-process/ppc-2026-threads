#pragma once

#include "romanova_v_linear_histogram_stretch/common/include/common.hpp"
#include "task/include/task.hpp"

namespace romanova_v_linear_histogram_stretch_threads {

class RomanovaVLinHistogramStretchSTL : public BaseTask {
 public:
  static constexpr ppc::task::TypeOfTask GetStaticTypeOfTask() {
    return ppc::task::TypeOfTask::kSTL;
  }
  explicit RomanovaVLinHistogramStretchSTL(const InType &in);

 private:
  bool ValidationImpl() override;
  bool PreProcessingImpl() override;
  bool RunImpl() override;
  bool PostProcessingImpl() override;
  void GetThreadRange(size_t thid, size_t total, size_t num_th, size_t &beg, size_t &en);
};

}  // namespace romanova_v_linear_histogram_stretch_threads
