#pragma once

#include <complex>
#include <vector>

#include "ermakov_a_spar_mat_mult/common/include/common.hpp"
#include "task/include/task.hpp"

namespace ermakov_a_spar_mat_mult {

class ErmakovASparMatMultOMP : public BaseTask {
 public:
  static constexpr ppc::task::TypeOfTask GetStaticTypeOfTask() {
    return ppc::task::TypeOfTask::kOMP;
  }
  explicit ErmakovASparMatMultOMP(const InType &in);

 private:
  bool ValidationImpl() override;
  bool PreProcessingImpl() override;
  bool RunImpl() override;
  bool PostProcessingImpl() override;

  static bool ValidateMatrix(const MatrixCRS &m);

  void ProcessRow(int i, std::vector<std::complex<double>> &row_vals, std::vector<int> &row_mark,
                  std::vector<int> &used_cols, std::vector<std::vector<std::complex<double>>> &row_values,
                  std::vector<std::vector<int>> &row_cols);

  void AccumulateRowProducts(int row_index, std::vector<std::complex<double>> &row_vals, std::vector<int> &row_mark,
                             std::vector<int> &used_cols);

  static void CollectRowValues(const std::vector<std::complex<double>> &row_vals, const std::vector<int> &used_cols,
                               std::vector<int> &cols, std::vector<std::complex<double>> &vals);

  static void SortUsedCols(std::vector<int> &cols);

  MatrixCRS a_;
  MatrixCRS b_;
  MatrixCRS c_;
};

}  // namespace ermakov_a_spar_mat_mult
