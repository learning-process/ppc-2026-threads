#pragma once

#include <vector>

#include "task/include/task.hpp"
#include "timur_a_cannon/common/include/common.hpp"

namespace timur_a_cannon {

class TimurACannonMatrixMultiplicationALL : public BaseTask {
 public:
  static constexpr ppc::task::TypeOfTask GetStaticTypeOfTask() {
    return ppc::task::TypeOfTask::kALL;
  }

  TimurACannonMatrixMultiplicationALL() = default;
  explicit TimurACannonMatrixMultiplicationALL(const InType &in);
  ~TimurACannonMatrixMultiplicationALL() override = default;

 private:
  bool ValidationImpl() override;
  bool PreProcessingImpl() override;
  bool RunImpl() override;
  bool PostProcessingImpl() override;

  static void BlockMultiplyAccumulate(const std::vector<std::vector<double>> &a,
                                      const std::vector<std::vector<double>> &b, std::vector<std::vector<double>> &c,
                                      int b_size);

  static std::vector<std::vector<double>> ComputeLocalResult(const std::vector<std::vector<double>> &src_a,
                                                             const std::vector<std::vector<double>> &src_b, int b_size,
                                                             int grid_sz, int block_row_start, int local_block_rows,
                                                             int n);
};

}  // namespace timur_a_cannon
