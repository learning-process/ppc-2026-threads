#pragma once

#include <cstdint>
#include <vector>

#include "buzulukski_d_gaus_gorizontal/common/include/common.hpp"
#include "task/include/task.hpp"

namespace buzulukski_d_gaus_gorizontal {

class BuzulukskiDGausGorizontalSTL : public BaseTask {
 public:
  static constexpr ppc::task::TypeOfTask GetStaticTypeOfTask() {
    return ppc::task::TypeOfTask::kSTL;
  }
  explicit BuzulukskiDGausGorizontalSTL(const InType &in);

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

  int width_ = 0;
  int height_ = 0;
  std::vector<uint8_t> input_image_;
  std::vector<uint8_t> output_image_;
};

}  // namespace buzulukski_d_gaus_gorizontal
