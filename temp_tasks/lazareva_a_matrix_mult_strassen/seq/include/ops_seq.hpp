#pragma once

#include <array>
#include <vector>

#include "lazareva_a_matrix_mult_strassen/common/include/common.hpp"
#include "task/include/task.hpp"

namespace lazareva_a_matrix_mult_strassen {

struct StrassenNode {
  std::vector<double> a;
  std::vector<double> b;
  int n{};
  int result_slot{};
  std::array<int, 7> child_slots{};
  bool expanded{};
};

class LazarevaATestTaskSEQ : public BaseTask {
 public:
  static constexpr ppc::task::TypeOfTask GetStaticTypeOfTask() {
    return ppc::task::TypeOfTask::kSEQ;
  }
  explicit LazarevaATestTaskSEQ(const InType &in);

 private:
  bool ValidationImpl() override;
  bool PreProcessingImpl() override;
  bool RunImpl() override;
  bool PostProcessingImpl() override;

  int n_{};
  int padded_n_{};
  std::vector<double> a_;
  std::vector<double> b_;
  std::vector<double> result_;

  static int NextPowerOfTwo(int n);
  static std::vector<double> PadMatrix(const std::vector<double> &m, int old_n, int new_n);
  static std::vector<double> UnpadMatrix(const std::vector<double> &m, int old_n, int new_n);
  static std::vector<double> Add(const std::vector<double> &a, const std::vector<double> &b, int n);
  static std::vector<double> Sub(const std::vector<double> &a, const std::vector<double> &b, int n);
  static void Split(const std::vector<double> &parent, int n, std::vector<double> &a11, std::vector<double> &a12,
                    std::vector<double> &a21, std::vector<double> &a22);
  static std::vector<double> Merge(const std::vector<double> &c11, const std::vector<double> &c12,
                                   const std::vector<double> &c21, const std::vector<double> &c22, int n);
  static std::vector<double> NaiveMult(const std::vector<double> &a, const std::vector<double> &b, int n);
  static std::vector<double> Strassen(const std::vector<double> &a, const std::vector<double> &b, int n);
};

}  // namespace lazareva_a_matrix_mult_strassen
