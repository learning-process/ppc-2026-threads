#pragma once

#include "task/include/task.hpp"
#include "volkov_a_sparse_mat_mul_ccs/common/include/common.hpp"

namespace volkov_a_sparse_mat_mul_ccs {

class VolkovASparseMatMulCcsStl : public BaseTask {
 public:
  static constexpr ppc::task::TypeOfTask GetStaticTypeOfTask() {
    return ppc::task::TypeOfTask::kSTL;
  }
  explicit VolkovASparseMatMulCcsStl(const InType &in);

 private:
  bool ValidationImpl() override;
  bool PreProcessingImpl() override;
  bool RunImpl() override;
  bool PostProcessingImpl() override;
};

}  // namespace volkov_a_sparse_mat_mul_ccs
