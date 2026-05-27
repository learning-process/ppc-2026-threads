#pragma once

#include "guseva_crs/common/include/common.hpp"
#include "task/include/task.hpp"

namespace guseva_crs {

class GusevaCRSMatMulStl : public BaseTask {
 public:
  static constexpr ppc::task::TypeOfTask GetStaticTypeOfTask() {
    return ppc::task::TypeOfTask::kSTL;
  }
  explicit GusevaCRSMatMulStl(const InType &in);

 private:
  bool ValidationImpl() override;
  bool PreProcessingImpl() override;
  bool RunImpl() override;
  bool PostProcessingImpl() override;
};

}  // namespace guseva_crs
