#pragma once

#include <vector>

#include "lobanov_d_multi_matrix_crs/common/include/common.hpp"
#include "task/include/task.hpp"

namespace lobanov_d_multi_matrix_crs {

class LobanovMultyMatrixALL : public BaseTask {
 public:
  static constexpr ppc::task::TypeOfTask GetStaticTypeOfTask() {
    return ppc::task::TypeOfTask::kALL;
  }
  explicit LobanovMultyMatrixALL(const InType &in);

 private:
  bool ValidationImpl() override;
  bool PreProcessingImpl() override;
  bool RunImpl() override;
  bool PostProcessingImpl() override;

  static void DistributeSparseMatrix(CompressedRowMatrix &mat, int root, int rows, int cols);
  static CompressedRowMatrix TransposeSparseMatrix(const CompressedRowMatrix &src);
  static void ComputeLocalProduct(const CompressedRowMatrix &a, const CompressedRowMatrix &b_tr, int start_row,
                                  int local_rows, std::vector<int> &row_nnz_counts, std::vector<double> &packed_vals,
                                  std::vector<int> &packed_cols);
  static void MergeLocalResults(int rank, int comm_size, int total_rows, int result_cols, int local_rows,
                                CompressedRowMatrix &result_mat, const std::vector<int> &row_nnz_counts,
                                const std::vector<double> &packed_vals, const std::vector<int> &packed_cols);
  static void SortIndices(std::vector<int> &vec);
};

}  // namespace lobanov_d_multi_matrix_crs
