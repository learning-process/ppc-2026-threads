#pragma once

#include "gasenin_l_djstra/common/include/common.hpp"
#include "task/include/task.hpp"

namespace gasenin_l_djstra {

class GaseninLDjstraALL : public BaseTask {
 public:
  static constexpr ppc::task::TypeOfTask GetStaticTypeOfTask() {
    return ppc::task::TypeOfTask::kALL;
  }
  explicit GaseninLDjstraALL(const InType &in);

 private:
  bool ValidationImpl() override;
  bool PreProcessingImpl() override;
  bool RunImpl() override;
  bool PostProcessingImpl() override;
};

}  // namespace gasenin_l_djstra
