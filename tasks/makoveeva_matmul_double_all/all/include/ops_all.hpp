#pragma once

#include <mpi.h>

#include <cstddef>
#include <vector>

#include "makoveeva_matmul_double_all/common/include/common.hpp"
#include "task/include/task.hpp"

namespace makoveeva_matmul_double_all {

class MatmulDoubleAllTask : public BaseTask {
 public:
  static constexpr ppc::task::TypeOfTask GetStaticTypeOfTask() {
    return ppc::task::TypeOfTask::kALL;
  }

  explicit MatmulDoubleAllTask(const InType &in);

  bool ValidationImpl() override;
  bool PreProcessingImpl() override;
  bool RunImpl() override;
  bool PostProcessingImpl() override;

  [[nodiscard]] const std::vector<double> &GetResult() const {
    return result_matrix_;
  }

  using BaseTask::GetInput;
  using BaseTask::GetOutput;

 private:
  void ParallelMultiply(size_t n, const std::vector<double> &a, const std::vector<double> &b, std::vector<double> &c);
  void SplitIntoBlocks(const std::vector<double> &src, std::vector<double> &dst, size_t n, size_t bs, int grid_size);
  void MergeFromBlocks(const std::vector<double> &src, std::vector<double> &dst, size_t n, size_t bs, int grid_size);
  void MultiplyBlockPair(const std::vector<double> &block_a, const std::vector<double> &block_b,
                         std::vector<double> &block_c, size_t bs);
  bool IsValidConfiguration(size_t n, int grid_size, int num_procs);
  void HandleFallback(int my_rank, size_t n, const std::vector<double> &a, const std::vector<double> &b,
                      std::vector<double> &c);
  void DistributeBlocks(int my_rank, const std::vector<double> &blocks_a, const std::vector<double> &blocks_b,
                        std::vector<double> &local_a, std::vector<double> &local_b, size_t block_sz);
  void ExecuteFoxIterations(int grid_dim, int row_id, int col_id, size_t bs, size_t block_sz, MPI_Comm row_comm,
                            std::vector<double> &local_a, std::vector<double> &local_b, std::vector<double> &local_c);
  void CollectResults(int my_rank, int num_procs, size_t n, size_t bs, size_t block_sz, int grid_dim,
                      const std::vector<double> &local_c, std::vector<double> &c);

  size_t matrix_size_ = 0;
  std::vector<double> matrix_a_;
  std::vector<double> matrix_b_;
  std::vector<double> result_matrix_;
};

}  // namespace makoveeva_matmul_double_all
