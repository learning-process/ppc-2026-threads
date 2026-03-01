#pragma once

#include <vector>

#include "task/include/task.hpp"
#include "zhurin_i_gauss_kernel_seq/common/include/common.hpp"

namespace zhurin_i_gauss_kernel_seq {

class ZhurinIGaussKernelSEQ : public BaseTask {
 public:
  static constexpr ppc::task::TypeOfTask GetStaticTypeOfTask() {
    return ppc::task::TypeOfTask::kSEQ;
  }

  explicit ZhurinIGaussKernelSEQ(const InType &in);

 private:
  bool ValidationImpl() override;
  bool PreProcessingImpl() override;
  bool RunImpl() override;
  bool PostProcessingImpl() override;

  static constexpr int kernel_[3][3] = {{1, 2, 1}, {2, 4, 2}, {1, 2, 1}};
  static constexpr int shift_ = 4;

  int width = 0;
  int height = 0;
  int numParts = 1;

  std::vector<std::vector<int>> image;
  std::vector<std::vector<int>> result;
  bool output_written = false;
};

}  // namespace zhurin_i_gauss_kernel_seq
