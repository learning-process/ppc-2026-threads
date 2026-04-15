#pragma once

#include "balchunayte_z_sobel/common/include/common.hpp"
#include "task/include/task.hpp"

namespace balchunayte_z_sobel {

class BalchunayteZSobelOpOMP : public BaseTask {
 public:
  static constexpr ppc::task::TypeOfTask GetStaticTypeOfTask() {
    return ppc::task::TypeOfTask::kOMP;
  }
  explicit BalchunayteZSobelOpOMP(const InType &input_image);

 private:
  bool ValidationImpl() override;
  bool PreProcessingImpl() override;
  bool RunImpl() override;
  bool PostProcessingImpl() override;
};

}  // namespace balchunayte_z_sobel
