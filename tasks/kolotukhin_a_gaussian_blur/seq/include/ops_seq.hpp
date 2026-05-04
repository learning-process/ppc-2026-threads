#pragma once

#include <cstdint>
#include <vector>

#include "kolotukhin_a_gaussian_blur/common/include/common.hpp"
#include "task/include/task.hpp"

namespace kolotukhin_a_gaussian_blur {

class KolotukhinAGaussinBlureSEQ : public BaseTask {
 public:
  static constexpr ppc::task::TypeOfTask GetStaticTypeOfTask() {
    return ppc::task::TypeOfTask::kSEQ;
  }
  explicit KolotukhinAGaussinBlureSEQ(const InType &in);

 private:
  bool ValidationImpl() override;
  bool PreProcessingImpl() override;
  bool RunImpl() override;
  bool PostProcessingImpl() override;

  [[nodiscard]] static std::uint8_t GetPixel(const std::vector<std::uint8_t> &pixel_data, int img_width, int img_height,
                                             int pos_x, int pos_y);
};

}  // namespace kolotukhin_a_gaussian_blur
