#pragma once

#include "korolev_k_matrix_mult/common/include/common.hpp"
#include "task/include/task.hpp"

namespace korolev_k_matrix_mult {

class KorolevKMatrixMultOMP : public BaseTask {
 public:
  static constexpr ppc::task::TypeOfTask GetStaticTypeOfTask() {
    return ppc::task::TypeOfTask::kOMP;
  }
  explicit KorolevKMatrixMultOMP(const InType &in);

 private:
  bool ValidationImpl() override;
  bool PreProcessingImpl() override;
  bool RunImpl() override;
  bool PostProcessingImpl() override;
};

}  // namespace korolev_k_matrix_mult
