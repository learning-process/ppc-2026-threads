#pragma once

#include <cstddef>
#include <vector>

#include "kulik_a_mat_mul_double_ccs/common/include/common.hpp"
#include "task/include/task.hpp"

namespace kulik_a_mat_mul_double_ccs {

class KulikAMatMulDoubleCcsOMP : public BaseTask {
 public:
  static constexpr ppc::task::TypeOfTask GetStaticTypeOfTask() {
    return ppc::task::TypeOfTask::kOMP;
  }
  explicit KulikAMatMulDoubleCcsOMP(const InType &in);

 private:
  bool ValidationImpl() override;
  bool PreProcessingImpl() override;
  bool RunImpl() override;
  bool PostProcessingImpl() override;
  static void ProcessColumn(size_t j, int tid, const CCS &a, const CCS &b,
                            std::vector<std::vector<double>> &thread_accum, std::vector<std::vector<bool>> &thread_nz,
                            std::vector<std::vector<size_t>> &thread_nnz_rows,
                            std::vector<std::vector<double>> &local_values,
                            std::vector<std::vector<size_t>> &local_rows);
};

}  // namespace kulik_a_mat_mul_double_ccs
