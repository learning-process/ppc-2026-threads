#pragma once

#include <vector>

#include "luzan_e_double_sparse_matrix_mult/common/include/common.hpp"
#include "task/include/task.hpp"

namespace luzan_e_double_sparse_matrix_mult {

class LuzanEDoubleSparseMatrixMultSTL : public BaseTask {
 public:
  static constexpr ppc::task::TypeOfTask GetStaticTypeOfTask() {
    return ppc::task::TypeOfTask::kSTL;
  }
  explicit LuzanEDoubleSparseMatrixMultSTL(const InType &in);

  static SparseMatrix CalcProdSTL(const SparseMatrix &a, const SparseMatrix &b);

  static void AssembleResult(SparseMatrix &c, unsigned cols, const std::vector<std::vector<double>> &values_per_col,
                             const std::vector<std::vector<unsigned>> &rows_per_col);
  static void ProcessColumn(const SparseMatrix &a, const SparseMatrix &b, unsigned b_col,
                            std::vector<std::vector<double>> &values_per_col,
                            std::vector<std::vector<unsigned>> &rows_per_col);

  static void CollectNonZeros(const std::vector<double> &tmp_col, unsigned b_col,
                              std::vector<std::vector<double>> &values_per_col,
                              std::vector<std::vector<unsigned>> &rows_per_col);

  static void AccumulateColumn(const SparseMatrix &a, const SparseMatrix &b, unsigned b_col,
                               std::vector<double> &tmp_col);

 private:
  bool ValidationImpl() override;
  bool PreProcessingImpl() override;
  bool RunImpl() override;
  bool PostProcessingImpl() override;
};

}  // namespace luzan_e_double_sparse_matrix_mult
