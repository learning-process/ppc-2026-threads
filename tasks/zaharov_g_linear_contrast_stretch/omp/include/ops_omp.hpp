#pragma once

#include "task/include/task.hpp"
#include "zaharov_g_linear_contrast_stretch/common/include/common.hpp"

namespace zaharov_g_linear_contrast_stretch {

class ZaharovGLinContrStrOMP : public BaseTask {
 public:
  static constexpr ppc::task::TypeOfTask GetStaticTypeOfTask() {
    return ppc::task::TypeOfTask::kOMP;
  }
  explicit ZaharovGLinContrStrOMP(const InType &in);

 private:
  bool ValidationImpl() override;
  bool PreProcessingImpl() override;
  bool RunImpl() override;
  bool PostProcessingImpl() override;
};

}  // namespace zaharov_g_linear_contrast_stretch
