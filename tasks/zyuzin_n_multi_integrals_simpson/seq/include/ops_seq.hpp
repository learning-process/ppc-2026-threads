#pragma once

#include "task/include/task.hpp"
#include "zyuzin_n_multi_integrals_simpson/common/include/common.hpp"

namespace zyuzin_n_multi_integrals_simpson {

class ZyuzinNSimpsonSEQ : public BaseTask {
 public:
  static constexpr ppc::task::TypeOfTask GetStaticTypeOfTask() {
    return ppc::task::TypeOfTask::kSEQ;
  }
  explicit ZyuzinNSimpsonSEQ(const InType &in);

 private:
  bool ValidationImpl() override;
  bool PreProcessingImpl() override;
  bool RunImpl() override;
  bool PostProcessingImpl() override;

  double ComputeSimpsonMultiDim();

  static double GetSimpsonWeight(int index, int n);

  double result_{0.0};
};

}  // namespace zyuzin_n_multi_integrals_simpson
