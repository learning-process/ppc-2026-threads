#pragma once

#include "sabutay_sparse_complex_ccs_mult/common/include/common.hpp"
#include "task/include/task.hpp"

namespace sabutay_sparse_complex_ccs_mult {

class SabutaySparseComplexCcsMultSEQ : public BaseTask {
 public:
  static constexpr ppc::task::TypeOfTask GetStaticTypeOfTask() {
    return ppc::task::TypeOfTask::kSEQ;
  }
  explicit SabutaySparseComplexCcsMultSEQ(const InType &in);

  static void SpMM(const CCS &a, const CCS &b, CCS &c);

 private:
  bool ValidationImpl() override;
  bool PreProcessingImpl() override;
  bool RunImpl() override;
  bool PostProcessingImpl() override;
};

}  // namespace sabutay_sparse_complex_ccs_mult
