#pragma once

#include "../../common/include/common.hpp"
#include "task/include/task.hpp"

namespace kaur_a_matrix_multiply {

class KaurAMatrixMultiplySTL : public BaseTask {
 public:
  static constexpr ppc::task::TypeOfTask GetStaticTypeOfTask() {
    return ppc::task::TypeOfTask::kSTL;
  }
  explicit KaurAMatrixMultiplySTL(const InType &in);

 private:
  bool ValidationImpl() override;
  bool PreProcessingImpl() override;
  bool RunImpl() override;
  bool PostProcessingImpl() override;
};

}  // namespace kaur_a_matrix_multiply
