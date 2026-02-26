#pragma once

#include "zhurin_i_gauss_kernel_seq/common/include/common.hpp"
#include "task/include/task.hpp"

namespace zhurin_i_gauss_kernel_seq {

class ZhurinIGaussKernelSEQ : public ppc::task::Task<InType, OutType> {
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

  int width = 0;
  int height = 0;
  int numParts = 0;                
  std::vector<std::vector<int>> image;
  std::vector<std::vector<int>> result;

  static constexpr int kernel[3][3] = {{1, 2, 1}, {2, 4, 2}, {1, 2, 1}};
  static constexpr int slip = 4;  
};

}  // namespace zhurin_i_gauss_kernel_seq