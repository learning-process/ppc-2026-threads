#pragma once

#include <mpi.h>

#include <vector>

#include "chyokotov_a_dense_matrix_mul_foxs_algorithm/common/include/common.hpp"
#include "task/include/task.hpp"

namespace chyokotov_a_dense_matrix_mul_foxs_algorithm {

class ChyokotovADenseMatMulFoxAlgorithmALL : public BaseTask {
 public:
  static constexpr ppc::task::TypeOfTask GetStaticTypeOfTask() {
    return ppc::task::TypeOfTask::kALL;
  }
  explicit ChyokotovADenseMatMulFoxAlgorithmALL(const InType &in);

 private:
  bool ValidationImpl() override;
  bool PreProcessingImpl() override;
  bool RunImpl() override;
  bool PostProcessingImpl() override;

  static int CalcPaddedSize(int n, int q);
  static void PadMatrix(const std::vector<double> &src, std::vector<double> &dst, int original_n, int padded_n);

  static void Multiply(const std::vector<double> &a_block, const std::vector<double> &b_block,
                       std::vector<double> &c_block, int block_size);

  static void DistributeData(MPI_Comm comm, int worker_rank, int worker_size, int q, int block_size,
                             const std::vector<double> &matrix_a_full, const std::vector<double> &matrix_b_full,
                             std::vector<double> &local_a, std::vector<double> &local_b);

  static void FoxAlgorithm(MPI_Comm comm, int worker_rank, int q, int block_size, std::vector<double> &local_a,
                           std::vector<double> &local_b, std::vector<double> &local_c);

  static void CollectResult(MPI_Comm comm, int worker_rank, int worker_size, int q, int block_size,
                            std::vector<double> &flat_result, const std::vector<double> &local_c);
};

}  // namespace chyokotov_a_dense_matrix_mul_foxs_algorithm
