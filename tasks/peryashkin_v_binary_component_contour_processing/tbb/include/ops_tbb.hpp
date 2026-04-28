#pragma once

#include "peryashkin_v_binary_component_contour_processing/common/include/common.hpp"
#include "task/include/task.hpp"

namespace peryashkin_v_binary_component_contour_processing {

class PeryashkinVBinaryComponentContourProcessingTBB : public BaseTask {
 public:
  static constexpr ppc::task::TypeOfTask GetStaticTypeOfTask() {
    return ppc::task::TypeOfTask::kTBB;
  }
  explicit PeryashkinVBinaryComponentContourProcessingTBB(const InType &in);

 private:
  bool ValidationImpl() override;
  bool PreProcessingImpl() override;
  bool RunImpl() override;
  bool PostProcessingImpl() override;

  OutType local_out_;
};

}  // namespace peryashkin_v_binary_component_contour_processing
