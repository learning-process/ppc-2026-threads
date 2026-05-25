#pragma once

#include <vector>

#include "lobanov_d_multi_matrix_crs/common/include/common.hpp"
#include "task/include/task.hpp"

namespace lobanov_d_multi_matrix_crs {

class LobanovMultyMatrixSTL : public BaseTask {
 public:
  static constexpr ppc::task::TypeOfTask GetStaticTypeOfTask() {
    return ppc::task::TypeOfTask::kSTL;
  }
  explicit LobanovMultyMatrixSTL(const InType &in);

 private:
  bool ValidationImpl() override;
  bool PreProcessingImpl() override;
  bool RunImpl() override;
  bool PostProcessingImpl() override;

  static void SortIndices(std::vector<int> &vec);
  static CompressedRowMatrix MultiplyMatrices(const CompressedRowMatrix &a, const CompressedRowMatrix &b);
};

}  // namespace lobanov_d_multi_matrix_crs
