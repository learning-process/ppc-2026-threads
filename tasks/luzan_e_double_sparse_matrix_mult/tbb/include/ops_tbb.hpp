#pragma once
#include <vector>

#include "luzan_e_double_sparse_matrix_mult/common/include/common.hpp"
#include "task/include/task.hpp"

namespace luzan_e_double_sparse_matrix_mult {

class LuzanEDoubleSparseMatrixMultTBB : public BaseTask {
 public:
  static constexpr ppc::task::TypeOfTask GetStaticTypeOfTask() {
    return ppc::task::TypeOfTask::kTBB;
  }
  explicit LuzanEDoubleSparseMatrixMultTBB(const InType &in);

  static SparseMatrix CalcProdTBB(const SparseMatrix &a, const SparseMatrix &b);

  static void MultiplyColumn(const SparseMatrix &a, const SparseMatrix &b, unsigned b_col,
                             std::vector<double> &tmp_col);

  static void CompressColumn(const std::vector<double> &tmp_col, std::vector<double> &values,
                             std::vector<unsigned> &rows);

  static void BuildResult(SparseMatrix &c, const std::vector<std::vector<double>> &values_per_col,
                          const std::vector<std::vector<unsigned>> &rows_per_col);

 private:
  bool ValidationImpl() override;
  bool PreProcessingImpl() override;
  bool RunImpl() override;
  bool PostProcessingImpl() override;
};

}  // namespace luzan_e_double_sparse_matrix_mult
