#pragma once

#include "tabalaev_a_matrix_mul_strassen/common/include/common.hpp"
#include "task/include/task.hpp"

namespace tabalaev_a_matrix_mul_strassen {

struct StrassenFrame {
  std::vector<double> mat_a;
  std::vector<double> mat_b;
  int n;
  int stage;
};

class TabalaevAMatrixMulStrassenSEQ : public BaseTask {
 public:
  static constexpr ppc::task::TypeOfTask GetStaticTypeOfTask() {
    return ppc::task::TypeOfTask::kSEQ;
  }
  explicit TabalaevAMatrixMulStrassenSEQ(const InType &in);

 private:
  bool ValidationImpl() override;
  bool PreProcessingImpl() override;
  bool RunImpl() override;
  bool PostProcessingImpl() override;

  std::vector<double> StrassenMultiply(const std::vector<double> &mat_a, const std::vector<double> &mat_b, int n);
  std::vector<double> Add(const std::vector<double> &mat_a, const std::vector<double> &mat_b);
  std::vector<double> Subtract(const std::vector<double> &mat_a, const std::vector<double> &mat_b);

  int a_rows_ = 0;
  int a_cols_b_rows_ = 0;
  int b_cols_ = 0;

  int padded_n_ = 0;

  std::vector<double> padded_a_;
  std::vector<double> padded_b_;
  std::vector<double> result_c_;
};

}  // namespace tabalaev_a_matrix_mul_strassen
