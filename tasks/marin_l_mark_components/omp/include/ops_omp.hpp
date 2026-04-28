#pragma once

#include <cstdint>
#include <vector>

#include "marin_l_mark_components/common/include/common.hpp"
#include "task/include/task.hpp"

namespace marin_l_mark_components {

class MarinLMarkComponentsOMP : public BaseTask {
 public:
  static constexpr ppc::task::TypeOfTask GetStaticTypeOfTask() {
    return ppc::task::TypeOfTask::kOMP;
  }

  explicit MarinLMarkComponentsOMP(const InType &in);

 private:
  bool ValidationImpl() override;
  bool PreProcessingImpl() override;
  bool RunImpl() override;
  bool PostProcessingImpl() override;

  static bool IsBinary(const Image &img);

  void FirstPassOMP();
  void MergeStripeBorders();
  void SecondPassOMP();
  void ConvertLabelsToOutput();

  std::vector<std::uint8_t> binary_;
  std::vector<int> labels_flat_;
  Labels labels_;
  std::vector<int> parent_;
  std::vector<int> stripe_offsets_;
  std::vector<int> stripe_used_counts_;
  std::vector<int> root_to_compact_;
  std::vector<int> root_generation_;
  int height_ = 0;
  int width_ = 0;
  int stripe_count_ = 1;
  int max_label_id_ = 0;
  int generation_id_ = 1;
};

}  // namespace marin_l_mark_components
