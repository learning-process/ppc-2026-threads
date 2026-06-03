#pragma once

#include "akhmetov_daniil_strassen_dense_double/common/include/common.hpp"
#include "task/include/task.hpp"

namespace akhmetov_daniil_strassen_dense_double {

class AkhmetovDStrassenDenseDoubleSTL : public BaseTask {
 public:
  static constexpr ppc::task::TypeOfTask GetStaticTypeOfTask() {
    return ppc::task::TypeOfTask::kSTL;
  }
  explicit AkhmetovDStrassenDenseDoubleSTL(const InType &in);

 protected:
  bool ValidationImpl() override;
  bool PreProcessingImpl() override;
  bool RunImpl() override;
  bool PostProcessingImpl() override;
};

}  // namespace akhmetov_daniil_strassen_dense_double
