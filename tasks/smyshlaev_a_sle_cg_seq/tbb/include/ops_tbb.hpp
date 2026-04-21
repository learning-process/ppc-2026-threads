#pragma once

#include <vector>

#include "smyshlaev_a_sle_cg_seq/common/include/common.hpp"
#include "task/include/task.hpp"

namespace smyshlaev_a_sle_cg_seq {

class SmyshlaevASleCgTaskTBB : public BaseTask {
 public:
  static constexpr ppc::task::TypeOfTask GetStaticTypeOfTask() {
    return ppc::task::TypeOfTask::kTBB;
  }
  explicit SmyshlaevASleCgTaskTBB(const InType &in);

 private:
  std::vector<double> flat_A_;
  bool RunSequential();
  bool RunParallel(int num_threads);
  bool ValidationImpl() override;
  bool PreProcessingImpl() override;
  bool RunImpl() override;
  bool PostProcessingImpl() override;
};

}  // namespace smyshlaev_a_sle_cg_seq
