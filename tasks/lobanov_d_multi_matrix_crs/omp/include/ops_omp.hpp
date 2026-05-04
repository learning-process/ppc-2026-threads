#pragma once

#include "lobanov_d_multi_matrix_crs/common/include/common.hpp"
#include "task/include/task.hpp"

namespace lobanov_d_multi_matrix_crs {

class LobanovMultyMatrixOMP : public BaseTask {
 public:
  static constexpr ppc::task::TypeOfTask GetStaticTypeOfTask() {
    return ppc::task::TypeOfTask::kOMP;
  }
  explicit LobanovMultyMatrixOMP(const InType &in);

 private:
  bool ValidationImpl() override;
  bool PreProcessingImpl() override;
  bool RunImpl() override;
  bool PostProcessingImpl() override;
};

}  // namespace lobanov_d_multi_matrix_crs
