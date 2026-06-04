#pragma once

#include <cstddef>

#include "task/include/task.hpp"
#include "zyuzin_n_multi_integrals_simpson/common/include/common.hpp"

namespace zyuzin_n_multi_integrals_simpson {

class ZyuzinNSimpsonALL : public BaseTask {
 public:
  static constexpr ppc::task::TypeOfTask GetStaticTypeOfTask() {
    return ppc::task::TypeOfTask::kALL;
  }
  explicit ZyuzinNSimpsonALL(const InType &in);

 private:
  bool ValidationImpl() override;
  bool PreProcessingImpl() override;
  bool RunImpl() override;
  bool PostProcessingImpl() override;

  static double GetSimpsonWeight(int index, int n);
  double ComputeLocalSimpsonSum(size_t begin, size_t end);

  double result_{0.0};
};

}  // namespace zyuzin_n_multi_integrals_simpson
