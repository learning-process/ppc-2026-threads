#pragma once

#include <complex>
#include <cstddef>
#include <vector>

#include "ermakov_a_spar_mat_mult/common/include/common.hpp"
#include "task/include/task.hpp"

namespace ermakov_a_spar_mat_mult {

class ErmakovASparMatMultSTL : public BaseTask {
 public:
  static constexpr ppc::task::TypeOfTask GetStaticTypeOfTask() {
    return ppc::task::TypeOfTask::kSTL;
  }
  explicit ErmakovASparMatMultSTL(const InType &in);

 private:
  struct RowData {
    std::vector<int> cols;
    std::vector<std::complex<double>> vals;
  };

  struct Workspace {
    std::vector<std::complex<double>> accum;
    std::vector<int> marks;
    std::vector<int> touched_cols;
  };

  bool ValidationImpl() override;
  bool PreProcessingImpl() override;
  bool RunImpl() override;
  bool PostProcessingImpl() override;

  static bool ValidateMatrix(const MatrixCRS &m);
  static int ResolveThreadCount(int row_count, std::size_t total_work);
  static std::vector<int> BuildRowCosts(const MatrixCRS &a, const MatrixCRS &b);
  static std::vector<int> BuildThreadBoundaries(const std::vector<int> &row_costs, int thread_count);
  static void ResetWorkspace(Workspace &workspace, int cols);
  void MultiplyRow(int row_index, Workspace &workspace, RowData &row_data) const;
  void FinalizeResult(const std::vector<RowData> &rows_data);

  MatrixCRS a_;
  MatrixCRS b_;
  MatrixCRS c_;
};

}  // namespace ermakov_a_spar_mat_mult
