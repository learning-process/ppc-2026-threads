#pragma once

#include "sinev_a_mult_matrix_fox_algorithm/common/include/common.hpp"
#include "task/include/task.hpp"

namespace sinev_a_mult_matrix_fox_algorithm {

class SinevAMultMatrixFoxAlgorithmOMP : public BaseTask {
 public:
  static constexpr ppc::task::TypeOfTask GetStaticTypeOfTask() {
    return ppc::task::TypeOfTask::kOMP;
  }
  explicit SinevAMultMatrixFoxAlgorithmOMP(const InType &in);

 private:
  bool ValidationImpl() override;
  bool PreProcessingImpl() override;
  bool RunImpl() override;
  bool PostProcessingImpl() override;

  void SimpleMultiply(size_t n, const std::vector<double>& A, 
                    const std::vector<double>& B, std::vector<double>& C);

};

}  // namespace sinev_a_mult_matrix_fox_algorithm
