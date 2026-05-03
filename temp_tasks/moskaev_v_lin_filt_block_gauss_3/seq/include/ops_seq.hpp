#pragma once

#include <cstdint>
#include <vector>

#include "moskaev_v_lin_filt_block_gauss_3/common/include/common.hpp"
#include "task/include/task.hpp"

namespace moskaev_v_lin_filt_block_gauss_3 {

class MoskaevVLinFiltBlockGauss3SEQ : public BaseTask {
 public:
  static constexpr ppc::task::TypeOfTask GetStaticTypeOfTask() {
    return ppc::task::TypeOfTask::kSEQ;
  }
  explicit MoskaevVLinFiltBlockGauss3SEQ(const InType &in);

 private:
  bool ValidationImpl() override;
  bool PreProcessingImpl() override;
  bool RunImpl() override;
  bool PostProcessingImpl() override;

  static void ApplyGaussianFilterToBlock(const std::vector<uint8_t> &input_block, std::vector<uint8_t> &output_block,
                                         int block_width, int block_height, int channels);

  ImageInfo image_info_;
  int block_size_{0};
};

}  // namespace moskaev_v_lin_filt_block_gauss_3
