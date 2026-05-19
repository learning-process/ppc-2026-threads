#pragma once

#include "ashihmin_d_mult_matr_crs/common/include/common.hpp"
#include "task/include/task.hpp"

namespace ashihmin_d_mult_matr_crs {

class AshihminDMultMatrCrsTBB : public BaseTask {
 public:
  static constexpr ppc::task::TypeOfTask GetStaticTypeOfTask() {
    return ppc::task::TypeOfTask::kTBB;
  }

  explicit AshihminDMultMatrCrsTBB(const InType &in);

 private:
  bool ValidationImpl() override;
  bool PreProcessingImpl() override;
  bool RunImpl() override;
  bool PostProcessingImpl() override;
};

}  // namespace ashihmin_d_mult_matr_crs
