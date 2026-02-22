#pragma once

#include "task/include/task.hpp"
#include "timur_a_Cannon/common/include/common.hpp"

namespace timur_a_Cannon {

class TimurACannonSEQ : public BaseTask {
 public:
  static constexpr ppc::task::TypeOfTask GetStaticTypeOfTask() {
    return ppc::task::TypeOfTask::kSEQ;
  }

  explicit TimurACannonSEQ(const InType &in);

 private:
  bool ValidationImpl() override;
  bool PreProcessingImpl() override;
  bool RunImpl() override;
  bool PostProcessingImpl() override;
};

}  // namespace timur_a_Cannon
