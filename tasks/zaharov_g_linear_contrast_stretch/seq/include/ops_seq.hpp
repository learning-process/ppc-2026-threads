#pragma once

#include "zaharov_g_linear_contrast_stretch/common/include/common.hpp"
#include "task/include/task.hpp"

namespace zaharov_g_linear_contrast_stretch {

class ZaharovGLinContrStrSEQ : public BaseTask {
 public:
  static constexpr ppc::task::TypeOfTask GetStaticTypeOfTask() {
    return ppc::task::TypeOfTask::kSEQ;
  }
  explicit ZaharovGLinContrStrSEQ(const InType &in);

 private:
  bool ValidationImpl() override;
  bool PreProcessingImpl() override;
  bool RunImpl() override;
  bool PostProcessingImpl() override;
};

}  // namespace zaharov_g_linear_contrast_stretch
