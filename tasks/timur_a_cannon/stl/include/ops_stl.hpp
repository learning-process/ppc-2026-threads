#pragma once

#include <vector>

#include "task/include/task.hpp"
#include "timur_a_cannon/common/include/common.hpp"

namespace timur_a_cannon {

class TimurACannonMatrixMultiplicationSTL : public BaseTask {
 public:
  static constexpr ppc::task::TypeOfTask GetStaticTypeOfTask() {
    return ppc::task::TypeOfTask::kSTL;
  }

  TimurACannonMatrixMultiplicationSTL() = default;
  explicit TimurACannonMatrixMultiplicationSTL(const InType &in);
  ~TimurACannonMatrixMultiplicationSTL() override = default;

 private:
  bool ValidationImpl() override;
  bool PreProcessingImpl() override;
  bool RunImpl() override;
  bool PostProcessingImpl() override;

  static void BlockMultiplyAccumulate(const std::vector<std::vector<double>> &a,
                                      const std::vector<std::vector<double>> &b, std::vector<std::vector<double>> &c,
                                      int b_size);
};

}  // namespace timur_a_cannon
