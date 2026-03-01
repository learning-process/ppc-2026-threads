#pragma once

#include "kopilov_d_vertical_gauss_filter/common/include/common.hpp"
#include "task/include/task.hpp"

namespace kopilov_d_vertical_gauss_filter {

class KopilovDVerticalGaussFilterSEQ : public BaseTask {
 public:
  static constexpr ppc::task::TypeOfTask GetStaticTypeOfTask() {
    return ppc::task::TypeOfTask::kSEQ;
  }
  explicit KopilovDVerticalGaussFilterSEQ(const InType &in);

 private:
  bool ValidationImpl() override;
  bool PreProcessingImpl() override;
  bool RunImpl() override;
  bool PostProcessingImpl() override;
};

}  // namespace kopilov_d_vertical_gauss_filter