#pragma once

#include "sinev_a_mult_matrix_fox_algorithm/common/include/common.hpp"
#include "task/include/task.hpp"

namespace sinev_a_mult_matrix_fox_algorithm_seq {

class SinevAMultMatrixFoxAlgorithmSEQ : public BaseTask {
 public:
  static constexpr ppc::task::TypeOfTask GetStaticTypeOfTask() {
    return ppc::task::TypeOfTask::kSEQ;
  }
  explicit SinevAMultMatrixFoxAlgorithmSEQ(const InType &in);

 private:
  bool ValidationImpl() override;
  bool PreProcessingImpl() override;
  bool RunImpl() override;
  bool PostProcessingImpl() override;
};

}  // namespace sinev_a_mult_matrix_fox_algorithm_seq
