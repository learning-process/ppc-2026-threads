#pragma once

#include <cstddef>
#include <vector>

#include "romanova_v_linear_histogram_stretch/common/include/common.hpp"
#include "task/include/task.hpp"

namespace romanova_v_linear_histogram_stretch_threads {

class RomanovaVLinHistogramStretchALL : public BaseTask {
 public:
  static constexpr ppc::task::TypeOfTask GetStaticTypeOfTask() {
    return ppc::task::TypeOfTask::kALL;
  }
  explicit RomanovaVLinHistogramStretchALL(const InType &in);

 private:
  bool ValidationImpl() override;
  bool PreProcessingImpl() override;
  bool RunImpl() override;
  bool PostProcessingImpl() override;

  InType local_data_;
  OutType local_out_;
  size_t local_size_{};
  std::vector<int> vector_displs_;
  std::vector<int> vector_counts_;
};

}  // namespace romanova_v_linear_histogram_stretch_threads
