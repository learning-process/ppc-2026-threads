#pragma once

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
  void ParallelMultiplyBlocks(const std::vector<double>& A, const std::vector<double>& B, std::vector<double>& C, int size_block);
  void DistributeData(int rank, int size, int q, int size_block, std::vector<double>& local_A, std::vector<double>& local_B);
  void CannonAlgorithm(int rank, int q, int size_block, std::vector<double>& local_A, std::vector<double>& local_B, std::vector<double>& local_C);
  void CollectAndBroadcast(int rank, int size, int q, int size_block, const std::vector<double>& local_C);
  
};

}  // namespace safronov_m_multiplication_matrix_blocksscheme_cannon
