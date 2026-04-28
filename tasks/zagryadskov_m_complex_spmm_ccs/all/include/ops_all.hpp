#pragma once

#include <complex>
#include <vector>

#include "task/include/task.hpp"
#include "zagryadskov_m_complex_spmm_ccs/common/include/common.hpp"

namespace zagryadskov_m_complex_spmm_ccs {

class ZagryadskovMComplexSpMMCCSALL : public BaseTask {
 public:
  static constexpr ppc::task::TypeOfTask GetStaticTypeOfTask() {
    return ppc::task::TypeOfTask::kALL;
  }
  explicit ZagryadskovMComplexSpMMCCSALL(const InType &in);

 private:
  inline static void SpMM(const CCS &a, const CCS &b, CCS &c);
  inline static void SpMMSymbolic(const CCS &a, const CCS &b, std::vector<int> &col_ptr, int jstart, int jend);
  inline static void SpMMNumeric(const CCS &a, const CCS &b, CCS &c, const std::complex<double> &zero, int jstart,
                                 int jend);
  inline static void SpMMKernel(const CCS &a, const CCS &b, CCS &c, const std::complex<double> &zero,
                                std::vector<int> &rows, std::vector<std::complex<double>> &acc,
                                std::vector<int> &marker, int j);
  static void BcastCCS(CCS &a, int rank);
  static void ScatterB(const CCS &b, CCS &b_local, const std::vector<int> &col_starts, int rank, int size);
  static void SendCCS(const CCS &m, int dest);
  static void RecvCCS(CCS &m, int src);
  static void GatherC(CCS &c, CCS &c_local, int rank, int size);
  bool ValidationImpl() override;
  bool PreProcessingImpl() override;
  bool RunImpl() override;
  bool PostProcessingImpl() override;
};

}  // namespace zagryadskov_m_complex_spmm_ccs
