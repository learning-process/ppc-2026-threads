#pragma once

#include <vector>

#include "fedoseev_linear_image_filtering_vertical/common/include/common.hpp"
#include "task/include/task.hpp"

namespace fedoseev_linear_image_filtering_vertical {

class LinearImageFilteringVerticalAll : public BaseTask {
 public:
  static constexpr ppc::task::TypeOfTask GetStaticTypeOfTask() {
    return ppc::task::TypeOfTask::kMPI;
  }
  explicit LinearImageFilteringVerticalAll(const InType &in);

 private:
  bool ValidationImpl() override;
  bool PreProcessingImpl() override;
  bool RunImpl() override;
  bool PostProcessingImpl() override;

  void DistributeData(int rank, int size, int &h);
  void LocalProcessing(int rank, int size, int num_threads);
  void GatherData(int rank, int size);

  std::vector<int> counts_;
  std::vector<int> displs_;
  std::vector<int> local_data_;
  std::vector<int> local_result_;
  int w_ = 0;
};

}  // namespace fedoseev_linear_image_filtering_vertical
