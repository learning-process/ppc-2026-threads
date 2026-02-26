#pragma once

#include "zhurin_i_gauss_kernel_seq/common/include/common.hpp"
#include "task/include/task.hpp"

namespace zhurin_i_test_task_threads {

class NesterovATestTaskTBB : public BaseTask {
 public:
  static constexpr ppc::task::TypeOfTask GetStaticTypeOfTask() {
    return ppc::task::TypeOfTask::kTBB;
  }
  explicit NesterovATestTaskTBB(const InType &in);

 private:
  bool ValidationImpl() override;
  bool PreProcessingImpl() override;
  bool RunImpl() override;
  bool PostProcessingImpl() override;
};

}  // namespace zhurin_i_test_task_threads
