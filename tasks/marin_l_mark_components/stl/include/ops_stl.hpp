#pragma once

#include <cstdint>
#include <vector>

#include "marin_l_mark_components/common/include/common.hpp"
#include "task/include/task.hpp"

namespace marin_l_mark_components {

class MarinLMarkComponentsSTL : public BaseTask {
 public:
  static constexpr ppc::task::TypeOfTask GetStaticTypeOfTask() {
    return ppc::task::TypeOfTask::kSTL;
  }
  explicit MarinLMarkComponentsSTL(const InType &in);

 private:
  bool ValidationImpl() override;
  bool PreProcessingImpl() override;
  bool RunImpl() override;
  bool PostProcessingImpl() override;

  static bool IsBinary(const Image &img);

  void FirstPassSTL();
  void MergeStripeBorders();
  void SecondPassSTL();
  void ConvertLabelsToOutput();

  std::vector<std::uint8_t> binary_;
  std::vector<int> labels_flat_;
  Labels labels_;
  std::vector<int> parent_;
  std::vector<int> stripe_bounds_;
  std::vector<int> stripe_base_label_;
  std::vector<int> stripe_used_label_end_;
  std::vector<int> root_to_compact_;
  int height_ = 0;
  int width_ = 0;
  int stripe_count_ = 1;
  int total_max_labels_ = 1;
};

}  // namespace marin_l_mark_components
