#pragma once

#include "task/include/task.hpp"
#include "tsyplakov_k_mul_double_crs_matrix/common/include/common.hpp"

namespace tsyplakov_k_mul_double_crs_matrix {

class TsyplakovKTestTaskOMP : public BaseTask {
 public:
  static constexpr ppc::task::TypeOfTask GetStaticTypeOfTask() {
    return ppc::task::TypeOfTask::kOMP;
  }

  explicit TsyplakovKTestTaskOMP(const InType &in);

 private:
  bool ValidationImpl() override;
  bool PreProcessingImpl() override;
  bool RunImpl() override;
  bool PostProcessingImpl() override;
};

}  // namespace tsyplakov_k_mul_double_crs_matrix
