#pragma once

#include <tbb/tbb.h>

#include <vector>

#include "ivanova_p_marking_components_on_binary_image/common/include/common.hpp"
#include "task/include/task.hpp"

namespace ivanova_p_marking_components_on_binary_image {

class IvanovaPMarkingComponentsOnBinaryImageALL : public BaseTask {
 public:
  static constexpr ppc::task::TypeOfTask GetStaticTypeOfTask() {
    return ppc::task::TypeOfTask::kALL;
  }
  explicit IvanovaPMarkingComponentsOnBinaryImageALL(const InType &in);

 private:
  bool ValidationImpl() override;
  bool PreProcessingImpl() override;
  bool RunImpl() override;
  bool PostProcessingImpl() override;

  Image input_image_;
  std::vector<int> labels_;
  std::vector<int> parent_;
  int current_label_ = 0;
  int width_ = 0;
  int height_ = 0;

  int FindRoot(int label);
  void UnionLabels(int label1, int label2);
  void ProcessPixel(int xx, int yy, int idx);
  void FirstPass();
  void SecondPass();

  void InitLabelsAll(int total_pixels);
  void MergeHorizontalPairsAll();
  void MergeVerticalPairsAll();
  void FinalizeRootsAll(int total_pixels);
  void NormalizeLabelsAll(int total_pixels);
};

}  // namespace ivanova_p_marking_components_on_binary_image
