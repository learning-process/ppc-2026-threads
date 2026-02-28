#pragma once

#include "ashihmin_d_mult_matr_crs/common/include/common.hpp"

namespace ashihmin_d_mult_matr_crs {

class AshihminDMultMatrCrsSEQ : public BaseTask {
 public:
  static constexpr ppc::task::TypeOfTask GetStaticTypeOfTask() {
    return ppc::task::TypeOfTask::kSEQ;
  }

  explicit AshihminDMultMatrCrsSEQ(const InType &input_matrices);

 private:
  bool ValidationImpl() override;
  bool PreProcessingImpl() override;
  bool RunImpl() override;
  bool PostProcessingImpl() override;
};

}  // namespace ashihmin_d_mult_matr_crs
