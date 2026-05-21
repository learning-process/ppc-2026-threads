#pragma once

#include <oneapi/tbb/parallel_invoke.h>

#include <cstddef>
#include <vector>

#include "lazareva_a_matrix_mult_strassen/common/include/common.hpp"
#include "task/include/task.hpp"

namespace lazareva_a_matrix_mult_strassen {

class LazarevaATestTaskALL : public BaseTask {
 public:
  static constexpr ppc::task::TypeOfTask GetStaticTypeOfTask() {
    return ppc::task::TypeOfTask::kALL;
  }
  explicit LazarevaATestTaskALL(const InType &in);

  static std::vector<double> Add(const std::vector<double> &a, const std::vector<double> &b, int n);
  static std::vector<double> Sub(const std::vector<double> &a, const std::vector<double> &b, int n);
  static void Split(const std::vector<double> &parent, int n, std::vector<double> &a11, std::vector<double> &a12,
                    std::vector<double> &a21, std::vector<double> &a22);
  static std::vector<double> Merge(const std::vector<double> &c11, const std::vector<double> &c12,
                                   const std::vector<double> &c21, const std::vector<double> &c22, int h);
  static std::vector<double> NaiveMult(const std::vector<double> &a, const std::vector<double> &b, int n);

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
  static std::vector<double> StrassenTBB(const std::vector<double> &a, const std::vector<double> &b, int n);
  static std::vector<double> StrassenALL(const std::vector<double> &a, const std::vector<double> &b, int n);

  static std::vector<double> StrassenMaster(const std::vector<double> &a, const std::vector<double> &b, int h,
                                            size_t h_sz, int size);
  static void StrassenWorker(int rank, int h, size_t h_sz, int size);
};

}  // namespace lazareva_a_matrix_mult_strassen
