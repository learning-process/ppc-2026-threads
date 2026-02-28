#pragma once

#include "chyokotov_a_dense_matrix_mul_foxs_algorithm/common/include/common.hpp"
#include "task/include/task.hpp"

namespace chyokotov_a_dense_matrix_mul_foxs_algorithm {

class ChyokotovADenseMatMulFoxAlgorithmSEQ : public BaseTask {
 public:
  static constexpr ppc::task::TypeOfTask GetStaticTypeOfTask() {
    return ppc::task::TypeOfTask::kSEQ;
  }
  explicit ChyokotovADenseMatMulFoxAlgorithmSEQ(const InType &in);

 private:
  bool ValidationImpl() override;
  bool PreProcessingImpl() override;
  bool RunImpl() override;
  bool PostProcessingImpl() override;

  int CalculateBlockSize(int n);
  int CountBlock(int n, int size);
  void Matmul(std::vector<double> &a, std::vector<double> &b, int n, int istart, int iend, int jstart, int jend,
              int kstart, int kend);
};

}  // namespace chyokotov_a_dense_matrix_mul_foxs_algorithm
