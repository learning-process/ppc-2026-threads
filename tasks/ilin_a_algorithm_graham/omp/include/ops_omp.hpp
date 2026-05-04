#pragma once

#include <vector>

#include "ilin_a_algorithm_graham/common/include/common.hpp"
#include "task/include/task.hpp"

namespace ilin_a_algorithm_graham {

class IlinAGrahamOMP : public BaseTask {
 public:
  static constexpr ppc::task::TypeOfTask GetStaticTypeOfTask() {
    return ppc::task::TypeOfTask::kOMP;
  }
  explicit IlinAGrahamOMP(const InType &in);

 private:
  bool ValidationImpl() override;
  bool PreProcessingImpl() override;
  bool RunImpl() override;
  bool PostProcessingImpl() override;

  std::vector<Point> points_;
  std::vector<Point> hull_;
};

}  // namespace ilin_a_algorithm_graham
