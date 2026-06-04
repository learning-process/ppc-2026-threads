#pragma once

#include <array>
#include <vector>

#include "boltenkov_s_gaussian_kernel/common/include/common.hpp"
#include "task/include/task.hpp"

namespace boltenkov_s_gaussian_kernel {

class BoltenkovSGaussianKernelALL : public BaseTask {
 public:
  static constexpr ppc::task::TypeOfTask GetStaticTypeOfTask() {
    return ppc::task::TypeOfTask::kALL;
  }

  explicit BoltenkovSGaussianKernelALL(const InType &in);

 private:
  std::array<std::array<int, 3>, 3> kernel_;
  int shift_ = 4;

  bool ValidationImpl() override;
  bool PreProcessingImpl() override;
  bool RunImpl() override;
  bool PostProcessingImpl() override;

  static bool IsValidSize(int n, int m);
  static void ComputeScatterParams(int n, int m, int size, int rows_per_proc, std::vector<int> &send_counts,
                                   std::vector<int> &displs);
  static void ComputeGatherDispls(int m, const std::vector<int> &gather_counts, std::vector<int> &recv_counts,
                                  std::vector<int> &recv_displs);
  static std::vector<int> ApplyGaussianFilterFlat(const std::vector<int> &local_halo_flat, int halo_rows,
                                                  int local_start_row, int local_rows, int m,
                                                  const std::array<std::array<int, 3>, 3> &kernel, int shift);
};

}  // namespace boltenkov_s_gaussian_kernel
