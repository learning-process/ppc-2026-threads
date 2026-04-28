#pragma once

#include <cstdint>
#include <vector>

#include "rysev_m_linear_filter_gauss_kernel/common/include/common.hpp"
#include "task/include/task.hpp"

namespace rysev_m_linear_filter_gauss_kernel {

class RysevMGaussFilterOMP : public BaseTask {
 public:
  static constexpr ppc::task::TypeOfTask GetStaticTypeOfTask() {
    return ppc::task::TypeOfTask::kOMP;
  }
  explicit RysevMGaussFilterOMP(const InType &in);

  [[nodiscard]] const std::vector<uint8_t> &GetInputImage() const {
    return input_image_;
  }
  [[nodiscard]] const std::vector<uint8_t> &GetOutputImage() const {
    return output_image_;
  }
  [[nodiscard]] int GetWidth() const {
    return width_;
  }
  [[nodiscard]] int GetHeight() const {
    return height_;
  }
  [[nodiscard]] int GetChannels() const {
    return channels_;
  }

 private:
  bool ValidationImpl() override;
  bool PreProcessingImpl() override;
  bool RunImpl() override;
  bool PostProcessingImpl() override;

  void ApplyKernelToChannel(int channel, int rows, int cols);

  std::vector<uint8_t> input_image_;
  std::vector<uint8_t> output_image_;
  int width_ = 0;
  int height_ = 0;
  int channels_ = 3;
};

}  // namespace rysev_m_linear_filter_gauss_kernel
