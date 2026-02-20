#pragma once

#include <complex>
#include <vector>

#include "task/include/task.hpp"
#include "zagryadskov_m_complex_spmm_ccs/common/include/common.hpp"

namespace zagryadskov_m_complex_spmm_ccs {

class ZagryadskovMComplexSpMMCCSOMP : public BaseTask {
 public:
  static constexpr ppc::task::TypeOfTask GetStaticTypeOfTask() {
    return ppc::task::TypeOfTask::kOMP;
  }
  explicit ZagryadskovMComplexSpMMCCSOMP(const InType &in);

 private:
  static void SpMM(const CCS &a, const CCS &b, CCS &c);
  static void SpMMkernel(const CCS &a, const CCS &b, const std::complex<double> &zero, double eps, int num_threads,
                         std::vector<std::vector<int>> &t_row_ind,
                         std::vector<std::vector<std::complex<double>>> &t_values,
                         std::vector<std::vector<int>> &t_col_ptr);
  bool ValidationImpl() override;
  bool PreProcessingImpl() override;
  bool RunImpl() override;
  bool PostProcessingImpl() override;
};

}  // namespace zagryadskov_m_complex_spmm_ccs
