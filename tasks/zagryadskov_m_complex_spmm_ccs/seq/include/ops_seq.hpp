#pragma once

#include "task/include/task.hpp"
#include "zagryadskov_m_complex_spmm_ccs/common/include/common.hpp"

namespace zagryadskov_m_complex_spmm_ccs {

class ZagryadskovMComplexSpMMCCSSEQ : public BaseTask {
 public:
  static constexpr ppc::task::TypeOfTask GetStaticTypeOfTask() {
    return ppc::task::TypeOfTask::kSEQ;
  }
  explicit ZagryadskovMComplexSpMMCCSSEQ(const InType &in);

  static void SpMM(const CCS &a, const CCS &b, CCS &c);

 private:
  bool ValidationImpl() override;
  bool PreProcessingImpl() override;
  bool RunImpl() override;
  bool PostProcessingImpl() override;
};

}  // namespace zagryadskov_m_complex_spmm_ccs
