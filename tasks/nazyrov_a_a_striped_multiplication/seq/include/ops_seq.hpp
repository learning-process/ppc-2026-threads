#pragma once

#include "../../common/include/common.hpp"

namespace nazyrov_a_a_striped_multiplication {

class StripedMultiplicationSEQ : public BaseTask {
 public:
  static constexpr ppc::task::TypeOfTask GetStaticTypeOfTask() {
    return ppc::task::TypeOfTask::kSEQ;
  }

  explicit StripedMultiplicationSEQ(const InType &in);
  ~StripedMultiplicationSEQ() override = default;

 protected:
  bool ValidationImpl() override;
  bool PreProcessingImpl() override;
  bool RunImpl() override;
  bool PostProcessingImpl() override;

 private:
  std::vector<double> A_, B_, C_;
  int n_;  // размер матрицы (n x n)
};

}  // namespace nazyrov_a_a_striped_multiplication