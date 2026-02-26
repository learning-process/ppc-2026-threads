#pragma once

#include "ovsyannikov_n_simpson_method/common/include/common.hpp"

namespace ovsyannikov_n_simpson_method {

class OvsyannikovNSimpsonMethodSEQ : public BaseTask {
 public:
  static constexpr ppc::task::TypeOfTask GetStaticTypeOfTask() {
    return ppc::task::TypeOfTask::kSEQ;
  }
  explicit OvsyannikovNSimpsonMethodSEQ(const InType &in);

 private:
  bool ValidationImpl() override;
  bool PreProcessingImpl() override;
  bool RunImpl() override;
  bool PostProcessingImpl() override;

  static double function(double x, double y);  // Функция для интегрирования
  InType params_{};
  OutType res_{};
};

}  // namespace ovsyannikov_n_simpson_method
