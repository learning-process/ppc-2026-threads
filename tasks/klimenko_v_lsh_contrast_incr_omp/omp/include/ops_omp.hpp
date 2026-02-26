#pragma once

#include "klimenko_v_lsh_contrast_incr_omp/common/include/common.hpp"
#include "task/include/task.hpp"

namespace klimenko_v_lsh_contrast_incr_omp {

class KlimenkoVLSHContrastIncrOMP : public BaseTask {
 public:
  static constexpr ppc::task::TypeOfTask GetStaticTypeOfTask() {
    return ppc::task::TypeOfTask::kOMP;
  }
  explicit KlimenkoVLSHContrastIncrOMP(const InType &in);

 private:
  bool ValidationImpl() override;
  bool PreProcessingImpl() override;
  bool RunImpl() override;
  bool PostProcessingImpl() override;
};

}  // namespace klimenko_v_lsh_contrast_incr_omp
