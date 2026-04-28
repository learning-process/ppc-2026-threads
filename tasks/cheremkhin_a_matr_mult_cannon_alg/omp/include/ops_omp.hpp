#pragma once

#include "cheremkhin_a_matr_mult_cannon_alg/common/include/common.hpp"
#include "task/include/task.hpp"

namespace cheremkhin_a_matr_mult_cannon_alg {

class CheremkhinAMatrMultCannonAlgOMP : public BaseTask {
 public:
  static constexpr ppc::task::TypeOfTask GetStaticTypeOfTask() {
    return ppc::task::TypeOfTask::kOMP;
  }
  explicit CheremkhinAMatrMultCannonAlgOMP(const InType &in);

 private:
  bool ValidationImpl() override;
  bool PreProcessingImpl() override;
  bool RunImpl() override;
  bool PostProcessingImpl() override;
};

}  // namespace cheremkhin_a_matr_mult_cannon_alg
