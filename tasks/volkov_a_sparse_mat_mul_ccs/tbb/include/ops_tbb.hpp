#pragma once

#include <vector>

#include "volkov_a_sparse_mat_mul_ccs/common/include/common.hpp"
#include "task/include/task.hpp"

namespace volkov_a_sparse_mat_mul_ccs {

class VolkovASparseMatMulCcsTbb : public ppc::core::Task {
 public:
  static constexpr ppc::core::TaskType GetStaticTypeOfTask() {
    return ppc::core::TaskType::TBB;
  }
  explicit VolkovASparseMatMulCcsTbb(const InType &in);

 private:
  bool ValidationImpl() override;
  bool PreProcessingImpl() override;
  bool RunImpl() override;
  bool PostProcessingImpl() override;
};

}  // namespace volkov_a_sparse_mat_mul_ccs