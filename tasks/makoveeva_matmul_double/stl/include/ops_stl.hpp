#pragma once

#include <cstddef>
#include <mutex>
#include <vector>

#include "makoveeva_matmul_double/stl/include/common.hpp"
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

 private:
  size_t n_ = 0;
  std::vector<double> A_;
  std::vector<double> B_;
  std::vector<double> C_;

  // Вспомогательный метод для простого умножения маленьких матриц
  bool RunSimpleMultiply();

  // Worker функция для обработки диапазона итераций в отдельном потоке
  void Worker(size_t start_step, size_t end_step, size_t grid_size, size_t block_size, std::mutex &write_mutex);
};

}  // namespace makoveeva_matmul_double_stl
