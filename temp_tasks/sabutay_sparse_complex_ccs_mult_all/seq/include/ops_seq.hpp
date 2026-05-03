#pragma once

#include "sabutay_sparse_complex_ccs_mult_all/common/include/common.hpp"
#include "task/include/task.hpp"

namespace sabutay_sparse_complex_ccs_mult_all {

class SabutaySparseComplexCcsMultFixSEQ : public BaseTask {
 public:
  static constexpr ppc::task::TypeOfTask GetStaticTypeOfTask() {
    return ppc::task::TypeOfTask::kSEQ;
  }
  explicit SabutaySparseComplexCcsMultFixSEQ(const InType &in);

 private:
  static void BuildProductMatrix(const CCS &left, const CCS &right, CCS &out);
  bool ValidationImpl() override;
  bool PreProcessingImpl() override;
  bool RunImpl() override;
  bool PostProcessingImpl() override;
};

}  // namespace sabutay_sparse_complex_ccs_mult_all
