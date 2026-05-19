#pragma once

#include <vector>

#include "kopilov_d_vertical_gauss_filter/common/include/common.hpp"
#include "task/include/task.hpp"

namespace kopilov_d_vertical_gauss_filter {

class KopilovDVerticalGaussFilterALL : public BaseTask {
 public:
  static constexpr ppc::task::TypeOfTask GetStaticTypeOfTask() {
    return ppc::task::TypeOfTask::kALL;
  }
  explicit KopilovDVerticalGaussFilterALL(const InType &in);

 private:
  bool ValidationImpl() override;
  bool PreProcessingImpl() override;
  bool RunImpl() override;
  bool PostProcessingImpl() override;
};

}  // namespace kopilov_d_vertical_gauss_filter
