#pragma once

#include "sannikov_i_integrals_rectangle_method/common/include/common.hpp"
#include "task/include/task.hpp"

namespace sannikov_i_integrals_rectangle_method {

class SannikovIIntegralsRectangleMethodOMP : public BaseTask {
 public:
  static constexpr ppc::task::TypeOfTask GetStaticTypeOfTask() {
    return ppc::task::TypeOfTask::kOMP;
  }
  explicit SannikovIIntegralsRectangleMethodOMP(const InType &in);

 private:
  bool ValidationImpl() override;
  bool PreProcessingImpl() override;
  bool RunImpl() override;
  bool PostProcessingImpl() override;
};

}  // namespace sannikov_i_integrals_rectangle_method
