#pragma once

#include <vector>

#include "kamaletdinov_r_bitwise_int/common/include/common.hpp"
#include "task/include/task.hpp"

namespace kamaletdinov_r_bitwise_int {

class KamaletdinovRBitwiseIntOMP : public BaseTask {
 public:
  static constexpr ppc::task::TypeOfTask GetStaticTypeOfTask() {
    return ppc::task::TypeOfTask::kOMP;
  }
  explicit KamaletdinovRBitwiseIntOMP(const InType &in);

 private:
  bool ValidationImpl() override;
  bool PreProcessingImpl() override;
  bool RunImpl() override;
  bool PostProcessingImpl() override;

  std::vector<int> data_;
};

}  // namespace kamaletdinov_r_bitwise_int
