#pragma once

#include "galkin_d_multidim_integrals_rectangles/common/include/common.hpp"
#include "task/include/task.hpp"

namespace galkin_d_multidim_integrals_rectangles {

class GalkinDMultidimIntegralsRectanglesSEQ : public BaseTask {
 public:
  static constexpr ppc::task::TypeOfTask GetStaticTypeOfTask() {
    return ppc::task::TypeOfTask::kSEQ;
  }
  explicit GalkinDMultidimIntegralsRectanglesSEQ(const InType &in);

 private:
  bool ValidationImpl() override;
  bool PreProcessingImpl() override;
  bool RunImpl() override;
  bool PostProcessingImpl() override;
};

}  // namespace galkin_d_multidim_integrals_rectangles
