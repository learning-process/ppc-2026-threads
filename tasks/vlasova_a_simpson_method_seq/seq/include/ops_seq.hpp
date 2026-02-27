#pragma once
#include <vector>
#include "vlasova_a_simpson_method_seq/common/include/common.hpp"
#include "task/include/task.hpp"

namespace vlasova_a_simpson_method_seq {

class VlasovaASimpsonMethodSEQ : public BaseTask {
 public:
  static constexpr ppc::task::TypeOfTask GetStaticTypeOfTask() {
    return ppc::task::TypeOfTask::kSEQ;
  }

  explicit VlasovaASimpsonMethodSEQ(const InType &in);

 private:
  bool ValidationImpl() override;
  bool PreProcessingImpl() override;
  bool RunImpl() override;
  bool PostProcessingImpl() override;
  
  // Рекурсивное вычисление интеграла
  double SimpsonRecursive(size_t dim, std::vector<double>& point);
  
  InType task_data_;
  double result_;
  std::vector<double> h_;  // шаги интегрирования
};

}  // namespace vlasova_a_simpson_method_seq