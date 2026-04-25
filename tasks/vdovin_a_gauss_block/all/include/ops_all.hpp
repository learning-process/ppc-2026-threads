#pragma once

#include <cstdint>
#include <vector>

#include "task/include/task.hpp"
#include "vdovin_a_gauss_block/common/include/common.hpp"

namespace vdovin_a_gauss_block {

class VdovinAGaussBlockALL : public BaseTask {
 public:
  static constexpr ppc::task::TypeOfTask GetStaticTypeOfTask() {
    return ppc::task::TypeOfTask::kALL;
  }
  explicit VdovinAGaussBlockALL(const InType &in);

  std::vector<uint8_t> &InputImage() {
    return input_image_;
  }
  std::vector<uint8_t> &OutputImage() {
    return output_image_;
  }

 private:
  bool ValidationImpl() override;
  bool PreProcessingImpl() override;
  bool RunImpl() override;
  bool PostProcessingImpl() override;

  static uint8_t ComputePixelChannel(const std::vector<uint8_t> &img, int width, int height, int py, int px, int ch);

  int width_ = 0;
  int height_ = 0;
  std::vector<uint8_t> input_image_;
  std::vector<uint8_t> output_image_;
};

}  // namespace vdovin_a_gauss_block
