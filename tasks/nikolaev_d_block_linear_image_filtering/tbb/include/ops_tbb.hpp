#pragma once

#include <cstdint>
#include <vector>

#include "nikolaev_d_block_linear_image_filtering/common/include/common.hpp"
#include "task/include/task.hpp"

namespace nikolaev_d_block_linear_image_filtering {

class NikolaevDBlockLinearImageFilteringTBB : public BaseTask {
 public:
  static constexpr ppc::task::TypeOfTask GetStaticTypeOfTask() {
    return ppc::task::TypeOfTask::kTBB;
  }
  explicit NikolaevDBlockLinearImageFilteringTBB(const InType &in);

 private:
  bool ValidationImpl() override;
  bool PreProcessingImpl() override;
  bool RunImpl() override;
  bool PostProcessingImpl() override;

  static std::uint8_t GetPixel(const std::vector<uint8_t> &data, int w, int h, int x, int y, int ch);
  static std::uint8_t ApplyKernel(const std::vector<uint8_t> &src, int w, int h, int nx, int ny, int ch);
};

}  // namespace nikolaev_d_block_linear_image_filtering
