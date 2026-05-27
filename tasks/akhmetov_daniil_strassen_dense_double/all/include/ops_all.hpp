#pragma once

#include "akhmetov_daniil_strassen_dense_double/common/include/common.hpp"
#include "task/include/task.hpp"

namespace akhmetov_daniil_strassen_dense_double {

class AkhmetovDStrassenDenseDoubleALL : public BaseTask {
 public:
  static constexpr ppc::task::TypeOfTask GetStaticTypeOfTask() {
    return ppc::task::TypeOfTask::kALL;
  }
  explicit AkhmetovDStrassenDenseDoubleALL(const InType &in);

 private:
  bool ValidationImpl() override;
  bool PreProcessingImpl() override;
  bool RunImpl() override;
  bool PostProcessingImpl() override;
};

}  // namespace akhmetov_daniil_strassen_dense_double
