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

  // Helper methods to reduce cognitive complexity
  void ProcessStripePixelAll(int xx, int yy, int idx, int start_row, std::vector<int> &local_parent, int &local_label);
  static int FindLocalRootAll(int label, const std::vector<int> &local_parent);
  void UnionLocalLabelsAll(int label1, int label2, std::vector<int> &local_parent);
  void InitializeLocalParent(std::vector<int> &local_parent, int max_labels);
  void ProcessStripe(int start_row, int end_row, std::vector<int> &local_parent, int &local_label);
  void MergeBoundaries(int num_threads, int rows_per_thread);
  void MergeLocalParents(const std::vector<std::vector<int>> &local_parents, const std::vector<int> &local_labels,
                         int num_threads, int total_pixels);
};

}  // namespace ivanova_p_marking_components_on_binary_image
