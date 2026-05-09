#pragma once

#include "romanov_a_gauss_block/common/include/common.hpp"
#include "task/include/task.hpp"

namespace romanov_a_gauss_block {

class RomanovAGaussBlockOMP : public BaseTask {
 public:
  static constexpr ppc::task::TypeOfTask GetStaticTypeOfTask() {
    return ppc::task::TypeOfTask::kOMP;
  }
  explicit RomanovAGaussBlockOMP(const InType &in);

 private:
  bool ValidationImpl() override;
  bool PreProcessingImpl() override;
  bool RunImpl() override;
  bool PostProcessingImpl() override;
};

}  // namespace romanov_a_gauss_block
