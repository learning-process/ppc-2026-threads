#pragma once
#include <vector>

#include "luzan_e_double_sparse_matrix_mult/common/include/common.hpp"
#include "task/include/task.hpp"

namespace luzan_e_double_sparse_matrix_mult {

class LuzanEDoubleSparseMatrixMultALL : public BaseTask {
 public:
  static constexpr ppc::task::TypeOfTask GetStaticTypeOfTask() {
    return ppc::task::TypeOfTask::kALL;
  }
  explicit LuzanEDoubleSparseMatrixMultALL(const InType &in);

  static SparseMatrix CalcProdMPIOMP(const SparseMatrix &a_in, const SparseMatrix &b_in);

  static void BroadcastMatrix(SparseMatrix &m, int root);

  static void BuildColDistribution(int b_cols, int nprocs, std::vector<int> &counts, std::vector<int> &displs);

  static void ComputeLocalCols(const SparseMatrix &a, const SparseMatrix &b, int col_start, int col_count,
                               std::vector<std::vector<double>> &values_per_col,
                               std::vector<std::vector<unsigned>> &rows_per_col);

  static void GatherFlatArrays(int rank, int nprocs, const std::vector<double> &local_vals,
                               const std::vector<unsigned> &local_rows, std::vector<double> &global_vals,
                               std::vector<unsigned> &global_rows);

  static void FlattenLocalCols(const std::vector<std::vector<double>> &values_per_col,
                               const std::vector<std::vector<unsigned>> &rows_per_col, std::vector<int> &col_nnz,
                               std::vector<double> &flat_vals, std::vector<unsigned> &flat_rows);

 private:
  bool ValidationImpl() override;
  bool PreProcessingImpl() override;
  bool RunImpl() override;
  bool PostProcessingImpl() override;
};

}  // namespace luzan_e_double_sparse_matrix_mult
