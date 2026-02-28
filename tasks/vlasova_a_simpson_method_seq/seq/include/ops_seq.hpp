#pragma once
#include <cstddef>
#include <vector>

#include "task/include/task.hpp"
#include "vlasova_a_simpson_method_seq/common/include/common.hpp"

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

  void NextIndex(std::vector<int> &Index);
  double GetWeight(const std::vector<int> &Index) const;
  std::vector<double> GetPoint(const std::vector<int> &Index) const;

  InType task_data_;
  double result_;
  std::vector<double> h_;        // шаги интегрирования
  std::vector<int> dimensions_;  // количество точек по каждому измерению (n[i] + 1
};

}  // namespace vlasova_a_simpson_method_seq
