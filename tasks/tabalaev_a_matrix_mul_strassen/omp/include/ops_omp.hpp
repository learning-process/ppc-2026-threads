#pragma once

#include <cstddef>
#include <vector>

#include "tabalaev_a_matrix_mul_strassen/common/include/common.hpp"
#include "task/include/task.hpp"

namespace tabalaev_a_matrix_mul_strassen {

struct StrassenFrameOMP {
  std::vector<double> mat_a;
  std::vector<double> mat_b;
  size_t n;
  int stage;
};

class TabalaevAMatrixMulStrassenOMP : public BaseTask {
 public:
  static constexpr ppc::task::TypeOfTask GetStaticTypeOfTask() {
    return ppc::task::TypeOfTask::kOMP;
  }
  explicit TabalaevAMatrixMulStrassenOMP(const InType &in);

 private:
  bool ValidationImpl() override;
  bool PreProcessingImpl() override;
  bool RunImpl() override;
  bool PostProcessingImpl() override;

  static std::vector<double> StrassenMultiply(const std::vector<double> &mat_a, const std::vector<double> &mat_b,
                                              size_t n);

  static std::vector<double> Add(const std::vector<double> &mat_a, const std::vector<double> &mat_b);

  static std::vector<double> Subtract(const std::vector<double> &mat_a, const std::vector<double> &mat_b);

  static std::vector<double> BaseMultiply(const std::vector<double> &mat_a, const std::vector<double> &mat_b, size_t n);

  static void SplitMatrix(const std::vector<double> &src, size_t n, std::vector<double> &c11, std::vector<double> &c12,
                          std::vector<double> &c21, std::vector<double> &c22);

  static std::vector<double> CombineMatrix(const std::vector<double> &c11, const std::vector<double> &c12,
                                           const std::vector<double> &c21, const std::vector<double> &c22, size_t n);

  size_t a_rows_ = 0;
  size_t a_cols_b_rows_ = 0;
  size_t b_cols_ = 0;

  size_t padded_n_ = 0;

  std::vector<double> padded_a_;
  std::vector<double> padded_b_;
  std::vector<double> result_c_;
};

}  // namespace tabalaev_a_matrix_mul_strassen
