#pragma once

#include "kulik_a_mat_mul_double_ccs/common/include/common.hpp"
#include "task/include/task.hpp"

namespace kulik_a_mat_mul_double_ccs {

class KulikAMatMulDoubleCcsSEQ : public BaseTask {
 public:
  static constexpr ppc::task::TypeOfTask GetStaticTypeOfTask() {
    return ppc::task::TypeOfTask::kSEQ;
  }
  explicit KulikAMatMulDoubleCcsSEQ(const InType &in);

 private:
  bool ValidationImpl() override;
  bool PreProcessingImpl() override;
  bool RunImpl() override;
  bool PostProcessingImpl() override;
};

}  // namespace kulik_a_mat_mul_double_ccs
