#pragma once

#include <cstdint>
#include <vector>

#include "kolotukhin_a_gaussian_blur/common/include/common.hpp"
#include "task/include/task.hpp"

namespace kolotukhin_a_gaussian_blur {

class KolotukhinAGaussinBlurALL : public BaseTask {
 public:
  static constexpr ppc::task::TypeOfTask GetStaticTypeOfTask() {
    return ppc::task::TypeOfTask::kALL;
  }
  explicit KolotukhinAGaussinBlurALL(const InType &in);

 private:
  bool ValidationImpl() override;
  bool PreProcessingImpl() override;
  bool RunImpl() override;
  bool PostProcessingImpl() override;

  int rank_ = -1;
  int proc_count_ = 0;
  std::vector<std::uint8_t> local_data_;
  int local_height_ = 0;
  int global_height_ = 0;
  int global_width_ = 0;
  void DistributeWork();
  void SendWorkData(int rows_per_process, int remainder);
  void ReceiveWorkData();
  void GatherResults();
  void GatherResultsWorker(int rows_per_process, int remainder);
  void GatherResultsRoot(int rows_per_process, int remainder);
};

}  // namespace kolotukhin_a_gaussian_blur
