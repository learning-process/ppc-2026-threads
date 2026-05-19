#pragma once

#include <array>
#include <vector>

#include "task/include/task.hpp"
#include "zhurin_i_gaus_kernel/common/include/common.hpp"

namespace zhurin_i_gaus_kernel {

class ZhurinIGausKernelTBB : public BaseTask {
 public:
  static constexpr ppc::task::TypeOfTask GetStaticTypeOfTask() {
    return ppc::task::TypeOfTask::kTBB;
  }

  explicit ZhurinIGausKernelTBB(const InType &in);

 private:
  bool ValidationImpl() override;
  bool PreProcessingImpl() override;
  bool RunImpl() override;
  bool PostProcessingImpl() override;

  static constexpr std::array<std::array<int, 3>, 3> kKernel = {{{1, 2, 1}, {2, 4, 2}, {1, 2, 1}}};
  static constexpr int kShift = 4;

  int width_ = 0;
  int height_ = 0;
  int num_parts_ = 1;

  std::vector<std::vector<int>> padded_;
  std::vector<std::vector<int>> result_;
};

}  // namespace zhurin_i_gaus_kernel
