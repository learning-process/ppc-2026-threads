#pragma once

#include <unordered_map>
#include <vector>

#include "ivanova_p_marking_components_on_binary_image/common/include/common.hpp"
#include "task/include/task.hpp"

namespace ivanova_p_marking_components_on_binary_image {

class IvanovaPMarkingComponentsOnBinaryImageSEQ : public BaseTask {
 public:
  static constexpr ppc::task::TypeOfTask GetStaticTypeOfTask() {
    return ppc::task::TypeOfTask::kSEQ;
  }
  explicit IvanovaPMarkingComponentsOnBinaryImageSEQ(const InType &in);

 private:
  bool ValidationImpl() override;
  bool PreProcessingImpl() override;
  bool RunImpl() override;
  bool PostProcessingImpl() override;

  Image input_image_;
  std::vector<int> labels_;
  std::unordered_map<int, int> parent_;
  int current_label_ = 0;
  int width_ = 0;
  int height_ = 0;

  void FirstPass();
  void SecondPass();
  int FindRoot(int label);
  void UnionLabels(int label1, int label2);
};

}  // namespace ivanova_p_marking_components_on_binary_image
