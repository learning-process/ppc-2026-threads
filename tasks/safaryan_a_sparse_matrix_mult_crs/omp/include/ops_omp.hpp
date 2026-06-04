#pragma once

#include "safaryan_a_sparse_matrix_mult_crs/common/include/common.hpp"
#include "task/include/task.hpp"

namespace safaryan_a_sparse_matrix_mult_crs {

class SafaryanATaskOMP : public BaseTask {
 public:
  static constexpr ppc::task::TypeOfTask GetStaticTypeOfTask() {
    return ppc::task::TypeOfTask::kOMP;
  }
  explicit SafaryanATaskOMP(const InType &in);

 private:
  bool ValidationImpl() override;
  bool PreProcessingImpl() override;
  bool RunImpl() override;
  bool PostProcessingImpl() override;

  static SparseMatrixCCS MultiplyMatrices(const SparseMatrixCCS &a, const SparseMatrixCCS &b);
  static bool IsMatrixValid(const SparseMatrixCCS &matrix);
};

}  // namespace safaryan_a_sparse_matrix_mult_crs
