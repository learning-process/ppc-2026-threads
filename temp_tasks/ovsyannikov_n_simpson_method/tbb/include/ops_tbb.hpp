#pragma once
#include "ovsyannikov_n_simpson_method/common/include/common.hpp"
#include "task/include/task.hpp"

namespace ovsyannikov_n_simpson_method {
class OvsyannikovNSimpsonMethodTBB : public BaseTask {
 public:
  static constexpr ppc::task::TypeOfTask GetStaticTypeOfTask() {
    return ppc::task::TypeOfTask::kTBB;
  }
  explicit OvsyannikovNSimpsonMethodTBB(const InType &in);

  bool ValidationImpl() override;
  bool PreProcessingImpl() override;
  bool RunImpl() override;
  bool PostProcessingImpl() override;

 private:
  static double Function(double x, double y);
  static double GetCoeff(int i, int n);
  InType params_ = {};
  OutType res_ = 0.0;
};
}  // namespace ovsyannikov_n_simpson_method
