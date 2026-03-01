#pragma once

#include "example_processes/common/include/common.hpp"
#include "task/include/task.hpp"

namespace chaschin_v_linear_image_filtration_seq {

class ChaschinVLinearFiltrationSEQ : public BaseTask {
 public:
  static constexpr ppc::task::TypeOfTask GetStaticTypeOfTask() {
    return ppc::task::TypeOfTask::kSEQ;
  }
  explicit ChaschinVLinearFiltrationSEQ(const InType &in);

 private:
  bool ValidationImpl() override;
  bool PreProcessingImpl() override;
  bool RunImpl() override;
  bool PostProcessingImpl() override;
};

}  // namespace chaschin_v_linear_image_filtration_seq
