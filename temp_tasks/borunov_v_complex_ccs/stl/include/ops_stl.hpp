#pragma once

#include "borunov_v_complex_ccs/common/include/common.hpp"
#include "task/include/task.hpp"

namespace borunov_v_complex_ccs {

class BorunovVComplexCcsSTL : public BaseTask {
 public:
  static constexpr ppc::task::TypeOfTask GetStaticTypeOfTask() {
    return ppc::task::TypeOfTask::kSTL;
  }
  explicit BorunovVComplexCcsSTL(const InType &in);

 private:
  bool ValidationImpl() override;
  bool PreProcessingImpl() override;
  bool RunImpl() override;
  bool PostProcessingImpl() override;
};

}  // namespace borunov_v_complex_ccs
