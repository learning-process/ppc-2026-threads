#pragma once
#include <vector>

#include "task/include/task.hpp"
#include "vlasova_a_simpson_method_seq/common/include/common.hpp"

namespace vlasova_a_simpson_method_seq {

class VlasovaASimpsonMethodSEQ : public BaseTask {
 public:
  static constexpr ppc::task::TypeOfTask GetStaticTypeOfTask() {
    return ppc::task::TypeOfTask::kSEQ;
  }

  explicit VlasovaASimpsonMethodSEQ(InType in);

 private:
  bool ValidationImpl() override;
  bool PreProcessingImpl() override;
  bool RunImpl() override;
  bool PostProcessingImpl() override;

  void Nextindex(std::vector<int> &index);
  void ComputeWeight(const std::vector<int> &index, double &weight) const;
  void ComputePoint(const std::vector<int> &index, std::vector<double> &point) const;

  InType task_data_;
  double result_ = 0.0;
  std::vector<double> h_;        // шаги интегрирования
  std::vector<int> dimensions_;  // количество точек по каждому измерению (n[i] + 1
};

}  // namespace vlasova_a_simpson_method_seq
