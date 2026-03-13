#pragma once

#include <cstdint>
#include <vector>

#include "bortsova_a_integrals_rectangle/common/include/common.hpp"
#include "task/include/task.hpp"

namespace bortsova_a_integrals_rectangle {

class BortsovaAIntegralsRectangleOMP : public BaseTask {
 public:
  static constexpr ppc::task::TypeOfTask GetStaticTypeOfTask() {
    return ppc::task::TypeOfTask::kOMP;
  }
  explicit BortsovaAIntegralsRectangleOMP(const InType &in);

 private:
  bool ValidationImpl() override;
  bool PreProcessingImpl() override;
  bool RunImpl() override;
  bool PostProcessingImpl() override;

  Function func_;
  int dims_ = 0;
  int num_steps_ = 0;
  int64_t total_points_ = 0;
  double volume_ = 0.0;
  std::vector<std::vector<double>> midpoints_;
};

}  // namespace bortsova_a_integrals_rectangle
