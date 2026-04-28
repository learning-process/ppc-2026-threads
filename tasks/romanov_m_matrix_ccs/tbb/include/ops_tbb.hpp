#pragma once

#include <cstddef>
#include <vector>

#include "romanov_m_matrix_ccs/common/include/common.hpp"
#include "task/include/task.hpp"

namespace romanov_m_matrix_ccs {

class RomanovMMatrixCCSTBB : public BaseTask {
 public:
  static constexpr ppc::task::TypeOfTask GetStaticTypeOfTask() {
    return ppc::task::TypeOfTask::kTBB;
  }
  explicit RomanovMMatrixCCSTBB(const InType &in);

 private:
  bool ValidationImpl() override;
  bool PreProcessingImpl() override;
  bool RunImpl() override;
  bool PostProcessingImpl() override;

  static void MultiplyColumn(size_t col_index, const MatrixCCS &a, const MatrixCCS &b, std::vector<double> &temp_v,
                             std::vector<size_t> &temp_r);
};

}  // namespace romanov_m_matrix_ccs
