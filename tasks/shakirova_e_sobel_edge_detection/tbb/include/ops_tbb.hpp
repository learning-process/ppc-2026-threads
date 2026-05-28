#pragma once

#include "shakirova_e_sobel_edge_detection/common/include/common.hpp"
#include "task/include/task.hpp"

namespace shakirova_e_sobel_edge_detection {

class ShakirovaESobelEdgeDetectionTBB : public BaseTask {
 public:
  static constexpr ppc::task::TypeOfTask GetStaticTypeOfTask() {
    return ppc::task::TypeOfTask::kTBB;
  }
  explicit ShakirovaESobelEdgeDetectionTBB(const InType &in);

 private:
  bool ValidationImpl() override;
  bool PreProcessingImpl() override;
  bool RunImpl() override;
  bool PostProcessingImpl() override;
};

}  // namespace shakirova_e_sobel_edge_detection
