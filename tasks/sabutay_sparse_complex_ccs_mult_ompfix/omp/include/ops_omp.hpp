#pragma once

#include "sabutay_sparse_complex_ccs_mult_ompfix/common/include/common.hpp"
#include "task/include/task.hpp"

namespace sabutay_sparse_complex_ccs_mult_ompfix {

class SabutaySparseComplexCcsMultOmpFix : public BaseTask {
 public:
  static constexpr ppc::task::TypeOfTask GetStaticTypeOfTask() {
    return ppc::task::TypeOfTask::kOMP;
  }
  explicit SabutaySparseComplexCcsMultOmpFix(const InType &in);
  static void BuildProductMatrix(const CCS &left, const CCS &right, CCS &out);

 private:
  bool ValidationImpl() override;
  bool PreProcessingImpl() override;
  bool RunImpl() override;
  bool PostProcessingImpl() override;
};

}  // namespace sabutay_sparse_complex_ccs_mult_ompfix
