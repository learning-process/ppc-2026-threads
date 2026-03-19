#pragma once

#include <cstddef>
#include <vector>

#include "makoveeva_matmul_double_omp/common/include/common.hpp"
#include "task/include/task.hpp"

namespace makoveeva_matmul_double_omp {

class MatmulDoubleOMPTask : public makoveeva_matmul_double_seq::BaseTask {
 public:
  static constexpr ppc::task::TypeOfTask GetStaticTypeOfTask() {
    return ppc::task::TypeOfTask::kOMP;
  }
  
  explicit MatmulDoubleOMPTask(const makoveeva_matmul_double_seq::InType& in);

  bool ValidationImpl() override;
  bool PreProcessingImpl() override;
  bool RunImpl() override;
  bool PostProcessingImpl() override;

  [[nodiscard]] const std::vector<double>& GetResult() const { return C_; }

 private:
  size_t n_;                    // размер матрицы
  std::vector<double> A_;        // матрица A
  std::vector<double> B_;        // матрица B
  std::vector<double> C_;        // результат
};

}  // namespace makoveeva_matmul_double_omp