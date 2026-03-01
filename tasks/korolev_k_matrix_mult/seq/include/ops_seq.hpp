#pragma once

#include "korolev_k_matrix_mult/common/include/common.hpp"
#include "task/include/task.hpp"

namespace korolev_k_matrix_mult {

class KorolevKMatrixMultSEQ : public BaseTask {
 public:
  static constexpr ppc::task::TypeOfTask GetStaticTypeOfTask() {
    return ppc::task::TypeOfTask::kSEQ;
  }
  explicit KorolevKMatrixMultSEQ(const InType &in);

 private:
  bool ValidationImpl() override;
  bool PreProcessingImpl() override;
  bool RunImpl() override;
  bool PostProcessingImpl() override;
};

}  // namespace korolev_k_matrix_mult
