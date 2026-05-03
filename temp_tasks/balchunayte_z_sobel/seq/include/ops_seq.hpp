#pragma once

#include "balchunayte_z_sobel/common/include/common.hpp"
#include "task/include/task.hpp"

namespace balchunayte_z_sobel {

class BalchunayteZSobelOpSEQ : public BaseTask {
 public:
  static constexpr ppc::task::TypeOfTask GetStaticTypeOfTask() {
    return ppc::task::TypeOfTask::kSEQ;
  }
  explicit BalchunayteZSobelOpSEQ(const InType &input_image);

 private:
  bool ValidationImpl() override;
  bool PreProcessingImpl() override;
  bool RunImpl() override;
  bool PostProcessingImpl() override;
};

}  // namespace balchunayte_z_sobel
