#pragma once

#include "klimenko_v_lsh_contrast_incr_seq/common/include/common.hpp"
#include "task/include/task.hpp"

namespace klimenko_v_lsh_contrast_incr_seq {

class KlimenkoVLSHContrastIncrSEQ : public BaseTask {
 public:
  static constexpr ppc::task::TypeOfTask GetStaticTypeOfTask() {
    return ppc::task::TypeOfTask::kSEQ;
  }
  explicit KlimenkoVLSHContrastIncrSEQ(const InType &in);

 private:
  bool ValidationImpl() override;
  bool PreProcessingImpl() override;
  bool RunImpl() override;
  bool PostProcessingImpl() override;
};

}  // namespace klimenko_v_lsh_contrast_incr_seq
