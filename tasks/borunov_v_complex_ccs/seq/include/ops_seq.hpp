#pragma once

#include "borunov_v_complex_ccs/common/include/common.hpp"
#include "task/include/task.hpp"

namespace borunov_v_complex_ccs_seq {

class BorunovVComplexCcsSEQ : public BaseTask {
 public:
  static constexpr ppc::task::TypeOfTask GetStaticTypeOfTask() {
    return ppc::task::TypeOfTask::kSEQ;
  }
  explicit BorunovVComplexCcsSEQ(const InType &in);

 private:
  bool ValidationImpl() override;
  bool PreProcessingImpl() override;
  bool RunImpl() override;
  bool PostProcessingImpl() override;
};

}  // namespace borunov_v_complex_ccs_seq
