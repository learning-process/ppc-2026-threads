#pragma once

#include "iskhakov_d_vertical_gauss_filter/common/include/common.hpp"
#include "task/include/task.hpp"

namespace iskhakov_d_vertical_gauss_filter {

class IskhakovDVerticalGaussFilterSEQ : public BaseTask {
 public:
  static constexpr ppc::task::TypeOfTask GetStaticTypeOfTask() {
    return ppc::task::TypeOfTask::kSEQ;
  }
  explicit IskhakovDVerticalGaussFilterSEQ(const InType &in);

 private:
  bool ValidationImpl() override;
  bool PreProcessingImpl() override;
  bool RunImpl() override;
  bool PostProcessingImpl() override;
};

}  // namespace iskhakov_d_vertical_gauss_filter
