#pragma once

#include "gasenin_l_djstra/common/include/common.hpp"
#include "task/include/task.hpp"

namespace gasenin_l_djstra {

class GaseninLDjstraSTL : public BaseTask {
 public:
  static constexpr ppc::task::TypeOfTask GetStaticTypeOfTask() {
    return ppc::task::TypeOfTask::kSTL;
  }
  explicit GaseninLDjstraSTL(const InType &in);

 private:
  bool ValidationImpl() override;
  bool PreProcessingImpl() override;
  bool RunImpl() override;
  bool PostProcessingImpl() override;
};

}  // namespace gasenin_l_djstra
