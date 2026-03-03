#pragma once

#include <complex>
#include <vector>

#include "ermakov_a_spar_mat_mult_seq/common/include/common.hpp"
#include "task/include/task.hpp"

namespace ermakov_a_spar_mat_mult_seq {

class ErmakovASparMatMultSEQ : public BaseTask {
 public:
  static constexpr ppc::task::TypeOfTask GetStaticTypeOfTask() {
    return ppc::task::TypeOfTask::kSEQ;
  }
  explicit ErmakovASparMatMultSEQ(const InType &in);

 private:
  bool ValidationImpl() override;
  bool PreProcessingImpl() override;
  bool RunImpl() override;
  bool PostProcessingImpl() override;

  static bool ValidateMatrix(const MatrixCRS &m);
  void ProcessRow(int i, std::vector<std::complex<double>> &row_vals, std::vector<int> &row_mark,
                  std::vector<int> &used_cols, int &nnz_so_far);

  MatrixCRS a_;
  MatrixCRS b_;
  MatrixCRS c_;
};

}  // namespace ermakov_a_spar_mat_mult_seq
