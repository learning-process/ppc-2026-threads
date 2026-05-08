#pragma once

#include <utility>
#include <vector>

#include "chyokotov_a_dense_matrix_mul_foxs_algorithm/common/include/common.hpp"
#include "task/include/task.hpp"

namespace chyokotov_a_dense_matrix_mul_foxs_algorithm {

class ChyokotovADenseMatMulFoxAlgorithmSTL : public BaseTask {
 public:
  static constexpr ppc::task::TypeOfTask GetStaticTypeOfTask() {
    return ppc::task::TypeOfTask::kSTL;
  }
  explicit ChyokotovADenseMatMulFoxAlgorithmSTL(const InType &in);

 private:
  bool ValidationImpl() override;
  bool PreProcessingImpl() override;
  bool RunImpl() override;
  bool PostProcessingImpl() override;

  static int CalculateBlockSize(int n);
  static int CountBlock(int n, int size);
  static std::vector<std::pair<int, int>> Blocks(int count_block);
  void Matmul(std::vector<double> &a, std::vector<double> &b, int n, int istart, int iend, int jstart, int jend,
              int kstart, int kend);
};

}  // namespace chyokotov_a_dense_matrix_mul_foxs_algorithm
