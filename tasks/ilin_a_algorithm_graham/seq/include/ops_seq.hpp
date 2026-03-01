#pragma once

#include "ilin_a_algorithm_graham/common/include/common.hpp"
#include "task/include/task.hpp"

namespace ilin_a_algorithm_graham {

class IlinAGrahamSEQ : public BaseTask {
 public:
  static constexpr ppc::task::TypeOfTask GetStaticTypeOfTask() {
    return ppc::task::TypeOfTask::kSEQ;
  }
  explicit IlinAGrahamSEQ(const InType &in);

 private:
  bool ValidationImpl() override;
  bool PreProcessingImpl() override;
  bool RunImpl() override;
  bool PostProcessingImpl() override;

  std::vector<Point> points_;
  std::vector<Point> hull_;
};

}  // namespace ilin_a_algorithm_graham