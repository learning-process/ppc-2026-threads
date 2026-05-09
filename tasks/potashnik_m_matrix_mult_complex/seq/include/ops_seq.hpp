#pragma once

#include "potashnik_m_matrix_mult_complex/common/include/common.hpp"
#include "task/include/task.hpp"

namespace potashnik_m_matrix_mult_complex {

class PotashnikMMatrixMultComplexSEQ : public BaseTask {
 public:
  static constexpr ppc::task::TypeOfTask GetStaticTypeOfTask() {
    return ppc::task::TypeOfTask::kSEQ;
  }
  explicit PotashnikMMatrixMultComplexSEQ(const InType &in);

 private:
  bool ValidationImpl() override;
  bool PreProcessingImpl() override;
  bool RunImpl() override;
  bool PostProcessingImpl() override;
};

}  // namespace potashnik_m_matrix_mult_complex
