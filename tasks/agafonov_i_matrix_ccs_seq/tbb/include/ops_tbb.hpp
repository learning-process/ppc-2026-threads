#pragma once

#include <cstddef>
#include <vector>

#include "agafonov_i_matrix_ccs_seq/common/include/common.hpp"
#include "task/include/task.hpp"

namespace agafonov_i_matrix_ccs_seq {

class AgafonovIMatrixCCSTBB : public BaseTask {
 public:
  static constexpr ppc::task::TypeOfTask GetStaticTypeOfTask() {
    return ppc::task::TypeOfTask::kTBB;
  }
  explicit AgafonovIMatrixCCSTBB(const InType &in);

  bool ValidationImpl() override;
  bool PreProcessingImpl() override;
  bool RunImpl() override;
  bool PostProcessingImpl() override;

 private:
  static void ProcessColumn(size_t j, const InType::first_type &a, const InType::second_type &b,
                            std::vector<double> &accumulator, std::vector<size_t> &active_rows,
                            std::vector<bool> &row_mask, std::vector<double> &local_v, std::vector<int> &local_r);
};

}  // namespace agafonov_i_matrix_ccs_seq
