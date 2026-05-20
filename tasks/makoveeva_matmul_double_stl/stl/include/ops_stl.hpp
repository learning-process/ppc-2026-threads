#pragma once

#include <cstddef>
#include <vector>

#include "makoveeva_matmul_double_stl/common/include/common.hpp"
#include "task/include/task.hpp"

namespace makoveeva_matmul_double_stl {

class MatmulDoubleSTLTask : public BaseTask {
 public:
  static constexpr ppc::task::TypeOfTask GetStaticTypeOfTask() {
    return ppc::task::TypeOfTask::kSTL;
  }

  explicit MatmulDoubleSTLTask(const InType &in);

  bool ValidationImpl() override;
  bool PreProcessingImpl() override;
  bool RunImpl() override;
  bool PostProcessingImpl() override;

  [[nodiscard]] const std::vector<double> &GetResult() const {
    return C_;
  }

  using BaseTask::GetInput;
  using BaseTask::GetOutput;
  // hii
 private:
  size_t n_{0};
  std::vector<double> A_;
  std::vector<double> B_;
  std::vector<double> C_;
};

}  // namespace makoveeva_matmul_double_stl
