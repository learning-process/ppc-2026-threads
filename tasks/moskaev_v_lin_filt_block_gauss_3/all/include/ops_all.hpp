#pragma once

#include "moskaev_v_lin_filt_block_gauss_3/common/include/common.hpp"
#include "task/include/task.hpp"

namespace moskaev_v_lin_filt_block_gauss_3 {

class MoskaevVLinFiltBlockGauss3ALL : public BaseTask {
 public:
  static constexpr ppc::task::TypeOfTask GetStaticTypeOfTask() {
    return ppc::task::TypeOfTask::kALL;
  }
  explicit MoskaevVLinFiltBlockGauss3ALL(const InType &in);

 private:
  bool ValidationImpl() override;
  bool PreProcessingImpl() override;
  bool RunImpl() override;
  bool PostProcessingImpl() override;

  int rank_{0};
  int num_procs_{1};
  int block_size_{64};
};

}  // namespace moskaev_v_lin_filt_block_gauss_3
