#pragma once

#include "rychkova_gauss/common/include/common.hpp"
#include "task/include/task.hpp"

namespace rychkova_gauss {

class RychkovaGaussTBB : public BaseTask {
 public:
  static constexpr ppc::task::TypeOfTask GetStaticTypeOfTask() {
    return ppc::task::TypeOfTask::kTBB;
  }
  explicit RychkovaGaussTBB(const InType &in);

 private:
  bool ValidationImpl() override;
  bool PreProcessingImpl() override;
  bool RunImpl() override;
  bool PostProcessingImpl() override;
};

}  // namespace rychkova_gauss
