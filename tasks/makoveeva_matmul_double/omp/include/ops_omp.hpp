#pragma once

#include <cstddef>
#include <vector>

#include "makoveeva_matmul_double/omp/include/common.hpp"
#include "task/include/task.hpp"

namespace makoveeva_matmul_double_omp {

class MatmulDoubleOMPTask : public BaseTask {
 public:
  static constexpr ppc::task::TypeOfTask GetStaticTypeOfTask() {
    return ppc::task::TypeOfTask::kOMP;
  }

  explicit MatmulDoubleOMPTask(const InType &in);

  bool ValidationImpl() override;
  bool PreProcessingImpl() override;
  bool RunImpl() override;
  bool PostProcessingImpl() override;

  [[nodiscard]] const std::vector<double> &GetResult() const {
    return c_;
  }

 private:
  size_t n_ = 0;
  std::vector<double> a_;
  std::vector<double> b_;
  std::vector<double> c_;

  bool RunSimpleMultiply();
};

}  // namespace makoveeva_matmul_double_omp
