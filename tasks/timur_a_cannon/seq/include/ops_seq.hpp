#pragma once
#include <vector>

#include "task/include/task.hpp"
#include "timur_a_cannon/common/include/common.hpp"

namespace timur_a_cannon {

class TimurACannonMatrixMultiplication : public BaseTask {
 public:
  static constexpr ppc::task::TypeOfTask GetStaticTypeOfTask() {
    return ppc::task::TypeOfTask::kSEQ;
  }
  explicit TimurACannonMatrixMultiplication(const InType &in);

 private:
  bool ValidationImpl() override;
  bool PreProcessingImpl() override;
  bool RunImpl() override;
  bool PostProcessingImpl() override;

  static void RotateBlocksA(std::vector<std::vector<std::vector<std::vector<double>>>> &blocks, int grid_sz);
  static void RotateBlocksB(std::vector<std::vector<std::vector<std::vector<double>>>> &blocks, int grid_sz);
  static void BlockMultiplyAccumulate(const std::vector<std::vector<double>> &a,
                                      const std::vector<std::vector<double>> &b, std::vector<std::vector<double>> &c,
                                      int b_size);

  static void DistributeData(const std::vector<std::vector<double>> &src_a,
                             const std::vector<std::vector<double>> &src_b,
                             std::vector<std::vector<std::vector<std::vector<double>>>> &bl_a,
                             std::vector<std::vector<std::vector<std::vector<double>>>> &bl_b, int b_size, int grid_sz);

  static void CollectResult(const std::vector<std::vector<std::vector<std::vector<double>>>> &bl_c,
                            std::vector<std::vector<double>> &res, int b_size, int grid_sz);
};

}  // namespace timur_a_cannon
