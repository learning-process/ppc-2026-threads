#pragma once

#include "sabutay_sparse_complex_ccs_mult_tbb/common/include/common.hpp"
#include "task/include/task.hpp"

namespace sabutay_sparse_complex_ccs_mult_tbb {

class SabutayASparseComplexCcsMultTBB : public BaseTask {
 public:
  static constexpr ppc::task::TypeOfTask GetStaticTypeOfTask() {
    return ppc::task::TypeOfTask::kTBB;
  }
  explicit SabutayASparseComplexCcsMultTBB(const InType &in);

 private:
  static SparseMatrixCCS MultiplyTbb(const SparseMatrixCCS &lhs, const SparseMatrixCCS &rhs);

  bool ValidationImpl() override;
  bool PreProcessingImpl() override;
  bool RunImpl() override;
  bool PostProcessingImpl() override;

  bool is_input_valid_ = false;
  SparseMatrixCCS lhs_;
  SparseMatrixCCS rhs_;
  SparseMatrixCCS result_;
};

}  // namespace sabutay_sparse_complex_ccs_mult_tbb
