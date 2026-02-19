#pragma once

#include "kolotukhin_a_gaussian_blur/common/include/common.hpp"
#include "task/include/task.hpp"

namespace kolotukhin_a_gaussian_blur {

class KolotukhinAGaussinBlureSEQ : public BaseTask {
 public:
  static constexpr ppc::task::TypeOfTask GetStaticTypeOfTask() {
    return ppc::task::TypeOfTask::kSEQ;
  }
  explicit KolotukhinAGaussinBlureSEQ(const InType &in);

 private:
  bool ValidationImpl() override;
  bool PreProcessingImpl() override;
  bool RunImpl() override;
  bool PostProcessingImpl() override;
};

}  // namespace kolotukhin_a_gaussian_blur
