#pragma once

#include "dilshodov_a_spmm_double_css/common/include/common.hpp"
#include "task/include/task.hpp"

namespace dilshodov_a_spmm_double_css {

class DilshodovASpmmDoubleCssSeq : public BaseTask {
 public:
  static constexpr ppc::task::TypeOfTask GetStaticTypeOfTask() {
    return ppc::task::TypeOfTask::kSEQ;
  }

  explicit DilshodovASpmmDoubleCssSeq(const InType &in);

 private:
  bool ValidationImpl() override;
  bool PreProcessingImpl() override;
  bool RunImpl() override;
  bool PostProcessingImpl() override;
};

}  // namespace dilshodov_a_spmm_double_css
