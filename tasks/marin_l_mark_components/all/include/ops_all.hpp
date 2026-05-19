#pragma once

#include <cstdint>
#include <vector>

#include "marin_l_mark_components/common/include/common.hpp"
#include "task/include/task.hpp"

namespace marin_l_mark_components {

class MarinLMarkComponentsALL : public BaseTask {
 public:
  static constexpr ppc::task::TypeOfTask GetStaticTypeOfTask() {
    return ppc::task::TypeOfTask::kALL;
  }
  explicit MarinLMarkComponentsALL(const InType &in);

 private:
  bool ValidationImpl() override;
  bool PreProcessingImpl() override;
  bool RunImpl() override;
  bool PostProcessingImpl() override;

  static bool IsBinary(const Image &img);
  void ConvertLabelsToOutput();

  std::vector<std::uint8_t> local_binary_flat_;
  std::vector<int> global_labels_flat_;
  Labels labels_out_;
  int height_ = 0;
  int width_ = 0;
  int rank_ = 0;
  int world_size_ = 1;
};

}  // namespace marin_l_mark_components
