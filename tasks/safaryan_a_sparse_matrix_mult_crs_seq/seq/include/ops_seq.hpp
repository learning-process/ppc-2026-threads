#pragma once

#include "safaryan_a_sparse_matrix_mult_crs_seq/common/include/common.hpp"
#include "task/include/task.hpp"

namespace safaryan_a_sparse_matrix_mult_crs_seq {

class SafaryanASparseMatrixMultCRSSeq : public BaseTask {
 public:
  static constexpr ppc::task::TypeOfTask GetStaticTypeOfTask() {
    return ppc::task::TypeOfTask::kSEQ;
  }

  explicit SafaryanASparseMatrixMultCRSSeq(const InType &in);

 private:
  bool ValidationImpl() override;
  bool PreProcessingImpl() override;
  bool RunImpl() override;
  bool PostProcessingImpl() override;

  static bool IsMatrixValid(const CRSMatrix &m);
};

}  // namespace safaryan_a_sparse_matrix_mult_crs_seq
