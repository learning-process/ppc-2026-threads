#pragma once

#include <cstddef>
#include <vector>

#include "makoveeva_matmul_double_all/common/include/common.hpp"
#include "task/include/task.hpp"

namespace makoveeva_matmul_double_all {

class MatmulDoubleAllTask : public BaseTask {
 public:
  static constexpr ppc::task::TypeOfTask GetStaticTypeOfTask() {
    return ppc::task::TypeOfTask::kALL;
  }

  explicit MatmulDoubleAllTask(const InType &in);

  bool ValidationImpl() override;
  bool PreProcessingImpl() override;
  bool RunImpl() override;
  bool PostProcessingImpl() override;

  [[nodiscard]] const std::vector<double> &GetResult() const {
    return result_matrix_;
  }

  using BaseTask::GetInput;
  using BaseTask::GetOutput;

 private:
  size_t matrix_size_ = 0;
  std::vector<double> matrix_a_;
  std::vector<double> matrix_b_;
  std::vector<double> result_matrix_;
};

}  // namespace makoveeva_matmul_double_all
