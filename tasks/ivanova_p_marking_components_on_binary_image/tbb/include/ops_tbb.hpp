#pragma once

#include <tbb/tbb.h>

#include <vector>

#include "ivanova_p_marking_components_on_binary_image/common/include/common.hpp"
#include "task/include/task.hpp"

namespace ivanova_p_marking_components_on_binary_image {

class IvanovaPMarkingComponentsOnBinaryImageTBB : public BaseTask {
 public:
  static constexpr ppc::task::TypeOfTask GetStaticTypeOfTask() {
    return ppc::task::TypeOfTask::kTBB;
  }
  explicit IvanovaPMarkingComponentsOnBinaryImageTBB(const InType &in);

 private:
  bool ValidationImpl() override;
  bool PreProcessingImpl() override;
  bool RunImpl() override;
  bool PostProcessingImpl() override;

  Image input_image_;
  std::vector<int> labels_;
  std::vector<int> parent_;
  int width_ = 0;
  int height_ = 0;
  int current_label_ = 0;

  int FindRoot(int label);
  void UnionLabels(int label1, int label2);
  void ProcessPixel(int xx, int yy, int idx);
  void FirstPass();
  void SecondPass();

  void InitLabelsTbb(int total_pixels);
  void MergeHorizontalPairsTbb();
  void MergeVerticalPairsTbb();
  void FinalizeRootsTbb(int total_pixels);
  void NormalizeLabelsTbb(int total_pixels);

  // Helper methods to reduce cognitive complexity
  void ProcessStripePixel(int xx, int yy, int idx, int start_row);
  int FindLocalRoot(int label);
};

}  // namespace ivanova_p_marking_components_on_binary_image
