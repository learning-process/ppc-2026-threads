// tasks/peryashkin_v_binary_component_contour_processing/omp/include/ops_omp.hpp
#pragma once

#include "peryashkin_v_binary_component_contour_processing/common/include/common.hpp"
#include "task/include/task.hpp"

namespace peryashkin_v_binary_component_contour_processing {

class PeryashkinVBinaryComponentContourProcessingOMP : public BaseTask {
 public:
  static constexpr ppc::task::TypeOfTask GetStaticTypeOfTask() {
    return ppc::task::TypeOfTask::kOMP;
  }
  explicit PeryashkinVBinaryComponentContourProcessingOMP(const InType &in);

 private:
  bool ValidationImpl() override;
  bool PreProcessingImpl() override;
  bool RunImpl() override;
  bool PostProcessingImpl() override;

  OutType local_out_;
};

}  // namespace peryashkin_v_binary_component_contour_processing
