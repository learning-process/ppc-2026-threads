#pragma once

#include <vector>

#include "ivanova_p_marking_components_on_binary_image/common/include/common.hpp"
#include "task/include/task.hpp"

namespace ivanova_p_marking_components_on_binary_image {

class IvanovaPMarkingComponentsOnBinaryImageSTL : public BaseTask {
 public:
  static constexpr ppc::task::TypeOfTask GetStaticTypeOfTask() {
    return ppc::task::TypeOfTask::kSTL;
  }
  explicit IvanovaPMarkingComponentsOnBinaryImageSTL(const InType &in);

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

  void InitLabelsStl(int total_pixels, int num_threads);
  void MergeHorizontalPairsStl(int num_threads);
  void MergeVerticalPairsStl(int num_threads);
  void FinalizeRootsStl(int total_pixels, int num_threads);
  void NormalizeLabelsStl(int total_pixels);

  // Helper methods to reduce cognitive complexity
  void ProcessStripePixelStl(int xx, int yy, int idx, int start_row, std::vector<int> &local_parent, int &local_label);
  int FindLocalRootStl(int label, const std::vector<int> &local_parent);
  void UnionLocalLabelsStl(int label1, int label2, std::vector<int> &local_parent);
};

}  // namespace ivanova_p_marking_components_on_binary_image
