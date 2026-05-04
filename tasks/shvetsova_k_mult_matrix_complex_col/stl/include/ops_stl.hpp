#pragma once

#include "../../common/include/common.hpp"
#include "task/include/task.hpp"

namespace shvetsova_k_mult_matrix_complex_col {

class ShvetsovaKMultMatrixComplexSTL : public BaseTask {
 public:
  static constexpr ppc::task::TypeOfTask GetStaticTypeOfTask() {
    return ppc::task::TypeOfTask::kSTL;
  }
  explicit ShvetsovaKMultMatrixComplexSTL(const InType &in);

 private:
  bool ValidationImpl() override;
  bool PreProcessingImpl() override;
  bool RunImpl() override;
  bool PostProcessingImpl() override;
};

}  // namespace shvetsova_k_mult_matrix_complex_col
