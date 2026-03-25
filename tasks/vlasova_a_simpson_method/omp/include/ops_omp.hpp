#pragma once

#include <vector>
#include <cstddef>

#include "task/include/task.hpp"
#include "vlasova_a_simpson_method/common/include/common.hpp"

namespace vlasova_a_simpson_method {

class VlasovaASimpsonMethodOMP : public BaseTask {
 public:
  static constexpr ppc::task::TypeOfTask GetStaticTypeOfTask() {
    return ppc::task::TypeOfTask::kOMP;
  }

  explicit VlasovaASimpsonMethodOMP(InType in);

 private:
  bool ValidationImpl() override;
  bool PreProcessingImpl() override;
  bool RunImpl() override;
  bool PostProcessingImpl() override;

  InType task_data_;
  double result_ = 0.0;
  std::vector<double> h_;                    // шаги интегрирования
  std::vector<int> dimensions_;              // количество точек по каждому измерению n[i] + 1
  size_t total_points_;                      // общее количество точек
  std::vector<std::vector<double>> weights_; // предвычисленные веса для каждого измерения
};

}  // namespace vlasova_a_simpson_method