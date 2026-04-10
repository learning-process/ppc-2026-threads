#pragma once

#include <vector>

#include "shakirova_e_sobel_edge_detection/common/include/common.hpp"
#include "task/include/task.hpp"

namespace shakirova_e_sobel_edge_detection {

class ShakirovaESobelEdgeDetectionOMP : public BaseTask {
 public:
  static constexpr ppc::task::TypeOfTask GetStaticTypeOfTask() {
    return ppc::task::TypeOfTask::kOMP;
  }
  explicit ShakirovaESobelEdgeDetectionOMP(const InType &in);

 private:
  bool ValidationImpl() override;
  bool PreProcessingImpl() override;
  bool RunImpl() override;
  bool PostProcessingImpl() override;

  int width_{0};
  int height_{0};
  std::vector<int> input_;
};

}  // namespace shakirova_e_sobel_edge_detection
