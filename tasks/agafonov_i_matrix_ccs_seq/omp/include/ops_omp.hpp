#pragma once

#include <cstddef>
#include <vector>

#include "agafonov_i_matrix_ccs_seq/common/include/common.hpp"
#include "task/include/task.hpp"

namespace agafonov_i_matrix_ccs_seq {

class AgafonovIMatrixCCSOMP : public BaseTask {
 public:
  static constexpr ppc::task::TypeOfTask GetStaticTypeOfTask() {
    return ppc::task::TypeOfTask::kOMP;
  }
  explicit AgafonovIMatrixCCSOMP(const InType &in);

 private:
  bool ValidationImpl() override;
  bool PreProcessingImpl() override;
  bool RunImpl() override;
  bool PostProcessingImpl() override;
  static void ProcessColumn(size_t j, const InType::first_type &a, const InType::second_type &b,
                            std::vector<double> &accumulator, std::vector<size_t> &active_rows,
                            std::vector<bool> &row_mask, std::vector<double> &local_v, std::vector<int> &local_r);
};

}  // namespace agafonov_i_matrix_ccs_seq
