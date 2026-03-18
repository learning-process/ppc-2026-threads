#pragma once
#include <vector>

#include "task/include/task.hpp"
#include "timur_a_cannon/common/include/common.hpp"

namespace timur_a_cannon {

class TimurACannonMatrixMultiplicationOMP : public ppc::core::BaseTask<InType, OutType> {
 public:
  explicit TimurACannonMatrixMultiplicationOMP(const InType &in);
  bool ValidationImpl() override;
  bool PreProcessingImpl() override;
  bool RunImpl() override;
  bool PostProcessingImpl() override;

  static TestType GetStaticTypeOfTask() {
    return std::make_tuple("timur_a_cannon_omp", 0, std::vector<std::vector<double>>{},
                           std::vector<std::vector<double>>{}, std::vector<std::vector<double>>{});
  }
};

}  // namespace timur_a_cannon
