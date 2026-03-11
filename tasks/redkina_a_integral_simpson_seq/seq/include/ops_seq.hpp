#pragma once

#include <functional>
#include <vector>

#include "redkina_a_integral_simpson_seq/common/include/common.hpp"
#include "task/include/task.hpp"

namespace redkina_a_integral_simpson_seq {

class RedkinaAIntegralSimpsonSEQ : public BaseTask {
 public:
  static constexpr ppc::task::TypeOfTask GetStaticTypeOfTask() {
    return ppc::task::TypeOfTask::kSEQ;
  }

  explicit RedkinaAIntegralSimpsonSEQ(const InType &in);

 private:
  bool ValidationImpl() override;
  bool PreProcessingImpl() override;
  bool RunImpl() override;
  bool PostProcessingImpl() override;

  std::function<double(const std::vector<double> &)> func_;
  std::vector<double> a_;
  std::vector<double> b_;
  std::vector<int> n_;
  double result_ = 0.0;
};

}  // namespace redkina_a_integral_simpson_seq
