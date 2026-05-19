#pragma once

#include <cstddef>
#include <cstdint>

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
  static void GetThreadRange(size_t thid, size_t total, size_t num_th, size_t &beg, size_t &en);
  static void FindMinMax(size_t begin, size_t end, uint8_t &curr_min, uint8_t &curr_max, const InType &in);
};

}  // namespace romanova_v_linear_histogram_stretch_threads
