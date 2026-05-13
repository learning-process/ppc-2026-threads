#pragma once

#include <cstdint>
#include <vector>

#include "moskaev_v_lin_filt_block_gauss_3/common/include/common.hpp"
#include "task/include/task.hpp"

namespace moskaev_v_lin_filt_block_gauss_3 {

class MoskaevVLinFiltBlockGauss3ALL : public BaseTask {
 public:
  static constexpr ppc::task::TypeOfTask GetStaticTypeOfTask() {
    return ppc::task::TypeOfTask::kALL;
  }
  explicit MoskaevVLinFiltBlockGauss3ALL(const InType &in);

  static void ApplyGaussianFilterToBlock(const std::vector<uint8_t> &input_block, std::vector<uint8_t> &output_block,
                                         int block_width, int block_height, int channels);

 private:
  bool ValidationImpl() override;
  bool PreProcessingImpl() override;
  bool RunImpl() override;
  bool PostProcessingImpl() override;

  ImageInfo image_info_;
  int block_size_{0};
  int rank_{0};
  int num_procs_{1};
};

}  // namespace moskaev_v_lin_filt_block_gauss_3
