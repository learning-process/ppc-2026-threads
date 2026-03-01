#pragma once

#include <cstddef>
#include <stack>
#include <vector>

#include "muhammadkhon_i_stressen_alg/common/include/common.hpp"
#include "task/include/task.hpp"

namespace muhammadkhon_i_stressen_alg {

struct StrassenFrame {
  std::vector<double> mat_a;
  std::vector<double> mat_b;
  size_t n;
  int stage;
};

class MuhammadkhonIStressenAlgSEQ : public BaseTask {
 public:
  static constexpr ppc::task::TypeOfTask GetStaticTypeOfTask() {
    return ppc::task::TypeOfTask::kSEQ;
  }
  explicit MuhammadkhonIStressenAlgSEQ(const InType &in);

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
  static void PushStrassenSubtasks(std::stack<StrassenFrame> &frames, const std::vector<double> &mat_a,
                                   const std::vector<double> &mat_b, size_t n);
  static std::vector<double> CombineStrassenResults(std::stack<std::vector<double>> &results, size_t n);

  size_t a_rows_ = 0;
  size_t a_cols_b_rows_ = 0;
  size_t b_cols_ = 0;

  size_t padded_n_ = 0;

  std::vector<double> padded_a_;
  std::vector<double> padded_b_;
  std::vector<double> result_c_;
};

}  // namespace muhammadkhon_i_stressen_alg
