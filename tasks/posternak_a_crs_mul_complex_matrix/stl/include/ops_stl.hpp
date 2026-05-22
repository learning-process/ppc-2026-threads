#pragma once

#include "posternak_a_crs_mul_complex_matrix/common/include/common.hpp"
#include "task/include/task.hpp"

namespace posternak_a_crs_mul_complex_matrix {

class PosternakACRSMulComplexMatrixSTL : public BaseTask {
 public:
  static constexpr ppc::task::TypeOfTask GetStaticTypeOfTask() {
    return ppc::task::TypeOfTask::kSTL;
  }
  explicit PosternakACRSMulComplexMatrixSTL(const InType &in);

 private:
  bool ValidationImpl() override;
  bool PreProcessingImpl() override;
  bool RunImpl() override;
  bool PostProcessingImpl() override;
};

}  // namespace posternak_a_crs_mul_complex_matrix
