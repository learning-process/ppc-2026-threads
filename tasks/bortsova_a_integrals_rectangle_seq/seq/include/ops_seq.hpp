#pragma once

#include <vector>

#include "bortsova_a_integrals_rectangle_seq/common/include/common.hpp"
#include "task/include/task.hpp"

namespace bortsova_a_integrals_rectangle_seq {

class BortsovaAIntegralsRectangleSEQ : public BaseTask {
 public:
  static constexpr ppc::task::TypeOfTask GetStaticTypeOfTask() {
    return ppc::task::TypeOfTask::kSEQ;
  }
  explicit BortsovaAIntegralsRectangleSEQ(const InType &in);

 private:
  bool ValidationImpl() override;
  bool PreProcessingImpl() override;
  bool RunImpl() override;
  bool PostProcessingImpl() override;

  Function func_;
  std::vector<double> lower_bounds_;
  std::vector<double> upper_bounds_;
  int num_steps_ = 0;
};

}  // namespace bortsova_a_integrals_rectangle_seq
