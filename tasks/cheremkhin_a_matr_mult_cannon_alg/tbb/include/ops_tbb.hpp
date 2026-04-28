#pragma once

#include "cheremkhin_a_matr_mult_cannon_alg/common/include/common.hpp"
#include "task/include/task.hpp"

namespace cheremkhin_a_matr_mult_cannon_alg {

class CheremkhinAMatrMultCannonAlgTBB : public BaseTask {
 public:
  static constexpr ppc::task::TypeOfTask GetStaticTypeOfTask() {
    return ppc::task::TypeOfTask::kTBB;
  }
  explicit CheremkhinAMatrMultCannonAlgTBB(const InType &in);

 private:
  bool ValidationImpl() override;
  bool PreProcessingImpl() override;
  bool RunImpl() override;
  bool PostProcessingImpl() override;
};

}  // namespace cheremkhin_a_matr_mult_cannon_alg
