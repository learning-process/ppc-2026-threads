#pragma once

#include <mpi.h>

#include <cstddef>
#include <vector>

#include "sinev_a_mult_matrix_fox_algorithm/common/include/common.hpp"
#include "task/include/task.hpp"

namespace sinev_a_mult_matrix_fox_algorithm {

class SinevAMultMatrixFoxAlgorithmALL : public BaseTask {
 public:
  static constexpr ppc::task::TypeOfTask GetStaticTypeOfTask() {
    return ppc::task::TypeOfTask::kALL;
  }

  explicit SinevAMultMatrixFoxAlgorithmALL(const InType &in);

 private:
  bool ValidationImpl() override;
  bool PreProcessingImpl() override;
  bool RunImpl() override;
  bool PostProcessingImpl() override;

  static void SimpleMultiply(size_t n, const std::vector<double> &a, const std::vector<double> &b,
                             std::vector<double> &c);

  static void DecomposeToBlocks(const std::vector<double> &src, std::vector<double> &dst, size_t n, size_t bs, int q);

  static void AssembleFromBlocks(const std::vector<double> &src, std::vector<double> &dst, size_t n, size_t bs, int q);

  static void LocalMatrixMultiply(const std::vector<double> &local_a, const std::vector<double> &local_b,
                                  std::vector<double> &local_c, size_t bs);

  static bool NeedFallback(size_t n, int q, int world_size);

  static void ExecuteFallback(int rank, size_t n, const std::vector<double> &a, const std::vector<double> &b,
                              std::vector<double> &c);

  static void ScatterBlocks(int rank, const std::vector<double> &blocks_a, const std::vector<double> &blocks_b,
                            std::vector<double> &local_a, std::vector<double> &local_b, size_t block_size);

  static void RunFoxStages(int q, int row, int col, size_t bs, size_t block_size, MPI_Comm row_comm,
                           std::vector<double> &local_a, std::vector<double> &local_b, std::vector<double> &local_c);

  static void GatherResult(int rank, int world_size, size_t n, size_t bs, size_t block_size, int q,
                           const std::vector<double> &local_c, std::vector<double> &c);
};

}  // namespace sinev_a_mult_matrix_fox_algorithm
