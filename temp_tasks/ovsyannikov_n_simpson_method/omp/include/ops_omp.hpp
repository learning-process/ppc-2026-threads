#pragma once

#include "ovsyannikov_n_simpson_method/common/include/common.hpp"
#include "task/include/task.hpp"

namespace ovsyannikov_n_simpson_method {

class OvsyannikovNSimpsonMethodOMP : public BaseTask {
 public:
  static constexpr ppc::task::TypeOfTask GetStaticTypeOfTask() {
    return ppc::task::TypeOfTask::kOMP;
  }
  explicit OvsyannikovNSimpsonMethodOMP(const InType &in);

 private:
  bool ValidationImpl() override;
  bool PreProcessingImpl() override;
  bool RunImpl() override;
  bool PostProcessingImpl() override;

  static double Function(double x, double y);  // Функция для интегрирования
  static double GetCoeff(int i, int n);
  InType params_{};
  OutType res_{};
};

}  // namespace ovsyannikov_n_simpson_method
