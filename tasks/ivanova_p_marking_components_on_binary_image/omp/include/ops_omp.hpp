#pragma once

#include <vector>

#include "ivanova_p_marking_components_on_binary_image/common/include/common.hpp"
#include "task/include/task.hpp"

namespace ivanova_p_marking_components_on_binary_image {

class IvanovaPMarkingComponentsOnBinaryImageOMP : public BaseTask {
 public:
  static constexpr ppc::task::TypeOfTask GetStaticTypeOfTask() {
    return ppc::task::TypeOfTask::kOMP;
  }
  explicit IvanovaPMarkingComponentsOnBinaryImageOMP(const InType &in);

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

  // Методы для DSU (на векторе)
  int FindRoot(int i);
  void UnionLabels(int i, int j);

  // OMP-хелперы для снижения cognitive complexity
  void InitLabelsOmp(int total_pixels, int n_threads);
  void MergeHorizontalPairsOmp(int n_threads);
  void MergeVerticalPairsOmp(int n_threads);
  void FinalizeRootsOmp(int total_pixels, int n_threads);
  void NormalizeLabelsOmp(int total_pixels, int n_threads);
  static void TouchFrameworkOmp();
};

}  // namespace ivanova_p_marking_components_on_binary_image
