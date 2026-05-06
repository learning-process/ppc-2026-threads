#pragma once

#include <cstddef>
#include <vector>

#include "kapanova_s_sparse_matrix_mult_ccs/common/include/common.hpp"
#include "task/include/task.hpp"

namespace kapanova_s_sparse_matrix_mult_ccs {

class KapanovaSSparseMatrixMultCCSOMP : public BaseTask {
 public:
  static constexpr ppc::task::TypeOfTask GetStaticTypeOfTask() {
    return ppc::task::TypeOfTask::kOMP;
  }
  explicit KapanovaSSparseMatrixMultCCSOMP(const InType &in);

 private:
  bool ValidationImpl() override;
  bool PreProcessingImpl() override;
  bool RunImpl() override;
  bool PostProcessingImpl() override;

  static void ProcessColumn(int j, const CCSMatrix &a, const CCSMatrix &b,
                            std::vector<std::vector<double>> &thread_accum,
                            std::vector<std::vector<bool>> &thread_row_mask,
                            std::vector<std::vector<size_t>> &thread_active_rows,
                            std::vector<std::vector<std::vector<size_t>>> &thread_col_rows,
                            std::vector<std::vector<std::vector<double>>> &thread_col_vals);

  static void ComputeColumnSizes(int num_threads, size_t cols,
                                 const std::vector<std::vector<std::vector<size_t>>> &thread_col_rows,
                                 std::vector<size_t> &col_sizes);

  static void MergeThreadResults(int num_threads, size_t cols, CCSMatrix &c,
                                 const std::vector<std::vector<std::vector<size_t>>> &thread_col_rows,
                                 const std::vector<std::vector<std::vector<double>>> &thread_col_vals);
};

}  // namespace kapanova_s_sparse_matrix_mult_ccs
