#pragma once

#include <vector>  // <array> удалён

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

  static constexpr int kKernel[3][3] = {{1, 2, 1}, {2, 4, 2}, {1, 2, 1}};
  static constexpr int kShift = 4;

  int width_ = 0;
  int height_ = 0;
  int num_parts_ = 1;

  std::vector<std::vector<int>> image_;
  std::vector<std::vector<int>> result_;
  bool output_written_ = false;
};

}  // namespace zhurin_i_gauss_kernel_seq
