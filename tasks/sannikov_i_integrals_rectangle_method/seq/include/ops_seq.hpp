#pragma once

#include "sannikov_i_integrals_rectangle_method/common/include/common.hpp"
#include "task/include/task.hpp"

namespace sannikov_i_integrals_rectangle_method {

class SannikovIIntegralsRectangleMethodSEQ : public BaseTask {
 public:
  static constexpr ppc::task::TypeOfTask GetStaticTypeOfTask() {
    return ppc::task::TypeOfTask::kSEQ;
  }
  explicit SannikovIIntegralsRectangleMethodSEQ(const InType &in);

 private:
  bool ValidationImpl() override;
  bool PreProcessingImpl() override;
  bool RunImpl() override;
  bool PostProcessingImpl() override;
};

}  // namespace sannikov_i_integrals_rectangle_method
