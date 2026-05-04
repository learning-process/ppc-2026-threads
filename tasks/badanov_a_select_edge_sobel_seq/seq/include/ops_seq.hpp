#pragma once

#include <array>
#include <cstdint>
#include <vector>

#include "badanov_a_select_edge_sobel_seq/common/include/common.hpp"
#include "task/include/task.hpp"

namespace badanov_a_select_edge_sobel_seq {

class BadanovASelectEdgeSobelSEQ : public BaseTask {
 public:
  static constexpr ppc::task::TypeOfTask GetStaticTypeOfTask() {
    return ppc::task::TypeOfTask::kSEQ;
  }
  explicit BadanovASelectEdgeSobelSEQ(const InType &in);

 private:
  bool ValidationImpl() override;
  bool PreProcessingImpl() override;
  bool RunImpl() override;
  bool PostProcessingImpl() override;

  void ApplySobelOperator(const std::vector<uint8_t> &input, std::vector<float> &magnitude, float &max_magnitude);
  void ComputeGradientAtPixel(const std::vector<uint8_t> &input, int row, int col, float &gradient_x,
                              float &gradient_y) const;
  void ApplyThreshold(const std::vector<float> &magnitude, float max_magnitude, std::vector<uint8_t> &output) const;

  static constexpr std::array<std::array<int, 3>, 3> kKernelX = {{{{-1, 0, 1}}, {{-2, 0, 2}}, {{-1, 0, 1}}}};

  static constexpr std::array<std::array<int, 3>, 3> kKernelY = {{{{-1, -2, -1}}, {{0, 0, 0}}, {{1, 2, 1}}}};

  int width_ = 0;
  int height_ = 0;
  int threshold_ = 50;
};

}  // namespace badanov_a_select_edge_sobel_seq
