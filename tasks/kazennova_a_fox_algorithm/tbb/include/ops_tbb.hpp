#pragma once

#include "kazennova_a_fox_algorithm/common/include/common.hpp"
#include "task/include/task.hpp"

namespace kazennova_a_fox_algorithm {

class KazennovaATestTaskTBB : public BaseTask {
 public:
  static constexpr ppc::task::TypeOfTask GetStaticTypeOfTask() {
    return ppc::task::TypeOfTask::kTBB;
  }

  explicit KazennovaATestTaskTBB(const InType &in);

 private:
  bool ValidationImpl() override;
  bool PreProcessingImpl() override;
  bool RunImpl() override;
  bool PostProcessingImpl() override;

  static constexpr int kBlockSize = 64;
};

}  // namespace kazennova_a_fox_algorithm
