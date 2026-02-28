#pragma once

#include "lobanov_d_multi_matrix_crs/common/include/common.hpp"
#include "task/include/task.hpp"

namespace lobanov_d_multi_matrix_crs {
class LobanovMultyMatrixSEQ : public BaseTask {
 public:
  static constexpr ppc::task::TypeOfTask GetStaticTypeOfTask() {
    return ppc::task::TypeOfTask::kSEQ;
  }
  explicit LobanovMultyMatrixSEQ(const InType &in);

 private:
  bool ValidationImpl() override;
  bool PreProcessingImpl() override;
  bool RunImpl() override;
  bool PostProcessingImpl() override;

  static void PerformMatrixMultiplication(const CompressedRowMatrix &first_matrix,
                                          const CompressedRowMatrix &second_matrix,
                                          CompressedRowMatrix &product_result);
};

}  // namespace lobanov_d_multi_matrix_crs