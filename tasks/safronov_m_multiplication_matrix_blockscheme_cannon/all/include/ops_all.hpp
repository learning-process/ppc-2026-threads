#pragma once

#include <mpi.h>

#include <vector>

#include "safronov_m_multiplication_matrix_blockscheme_cannon/common/include/common.hpp"
#include "task/include/task.hpp"

namespace safronov_m_multiplication_matrix_blocksscheme_cannon {

class SafronovMMultiplicationMatrixBlockSchemeCannonALL : public BaseTask {
 public:
  static constexpr ppc::task::TypeOfTask GetStaticTypeOfTask() {
    return ppc::task::TypeOfTask::kALL;
  }

  explicit SafronovMMultiplicationMatrixBlockSchemeCannonALL(const InType &in);

 private:
  bool ValidationImpl() override;
  bool PreProcessingImpl() override;
  bool RunImpl() override;
  bool PostProcessingImpl() override;

  static int CalcPaddedSize(int n, int q);

  static void PadMatrix(const std::vector<std::vector<double>> &src, std::vector<std::vector<double>> &dst,
                        int padded_n);

  static void ParallelMultiplyBlocks(const std::vector<double> &a, const std::vector<double> &b, std::vector<double> &c,
                                     int block_size);

  static void DistributeData(MPI_Comm comm, int worker_rank, int worker_size, int q, int block_size,
                             const std::vector<std::vector<double>> &matrix_a_full,
                             const std::vector<std::vector<double>> &matrix_b_full, std::vector<double> &local_a,
                             std::vector<double> &local_b);

  static void CannonAlgorithm(MPI_Comm comm, int worker_rank, int q, int block_size, std::vector<double> &local_a,
                              std::vector<double> &local_b, std::vector<double> &local_c);

  static void CollectResult(MPI_Comm comm, int worker_rank, int worker_size, int q, int block_size,
                            std::vector<double> &flat_result, const std::vector<double> &local_c);
  static void FillResultFromBuffer(std::vector<double> &flat_result, const std::vector<double> &buffer, int row,
                                   int col, int block_size, int padded_n);  //
};

}  // namespace safronov_m_multiplication_matrix_blocksscheme_cannon
