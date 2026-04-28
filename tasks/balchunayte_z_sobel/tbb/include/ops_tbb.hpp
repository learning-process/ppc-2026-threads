#pragma once

#include "balchunayte_z_sobel/common/include/common.hpp"
#include "task/include/task.hpp"

namespace balchunayte_z_sobel {

class BalchunayteZSobelOpTBB : public BaseTask {
 public:
  static constexpr ppc::task::TypeOfTask GetStaticTypeOfTask() {
    return ppc::task::TypeOfTask::kTBB;
  }
  explicit BalchunayteZSobelOpTBB(const InType &input_image);

 private:
  bool ValidationImpl() override;
  bool PreProcessingImpl() override;
  bool RunImpl() override;
  bool PostProcessingImpl() override;
};

}  // namespace balchunayte_z_sobel
