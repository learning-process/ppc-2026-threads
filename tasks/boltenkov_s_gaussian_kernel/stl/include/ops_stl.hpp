#pragma once

#include <cstddef>
#include <vector>

#include "boltenkov_s_gaussian_kernel/common/include/common.hpp"
#include "task/include/task.hpp"

namespace boltenkov_s_gaussian_kernel {

class BoltenkovSGaussianKernelSTL : public BaseTask {
 public:
  static constexpr ppc::task::TypeOfTask GetStaticTypeOfTask() {
    return ppc::task::TypeOfTask::kSTL;
  }
  explicit BoltenkovSGaussianKernelSTL(const InType &in);

 private:
  std::vector<std::vector<int>> kernel_;
  int shift_ = 4;
  bool ValidationImpl() override;
  bool PreProcessingImpl() override;
  bool RunImpl() override;
  bool PostProcessingImpl() override;
  static void CopyRowsToPaddedBuffer(const std::vector<std::vector<int>> &data, std::vector<std::vector<int>> &tmp_data,
                                     std::size_t n, unsigned int num_threads);
  void ApplyGaussianKernel(const std::vector<std::vector<int>> &tmp_data, std::vector<std::vector<int>> &res,
                           std::size_t n, std::size_t m, unsigned int num_threads);
};

}  // namespace boltenkov_s_gaussian_kernel
