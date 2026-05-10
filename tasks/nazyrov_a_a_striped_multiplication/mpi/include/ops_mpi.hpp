#pragma once

#include "../../common/include/common.hpp"

namespace nazyrov_a_a_striped_multiplication {

class StripedMultiplicationMPI : public BaseTask {
 public:
  static constexpr ppc::task::TypeOfTask GetStaticTypeOfTask() {
    return ppc::task::TypeOfTask::kMPI;
  }

  explicit StripedMultiplicationMPI(const InType &in);
  ~StripedMultiplicationMPI() override = default;

 protected:
  bool ValidationImpl() override;
  bool PreProcessingImpl() override;
  bool RunImpl() override;
  bool PostProcessingImpl() override;

 private:
  std::vector<double> A_, B_, local_C_;
  int n_, rank_, size_;
  int rows_per_proc_;
};

}  // namespace nazyrov_a_a_striped_multiplication