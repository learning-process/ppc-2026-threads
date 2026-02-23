#pragma once

#include <vector>

#include "kamaletdinov_r_bitwise_int_seq/common/include/common.hpp"
#include "task/include/task.hpp"

namespace kamaletdinov_r_bitwise_int_seq {

void BitwiseSort(std::vector<int> &arr);

class KamaletdinovRBitwiseIntSEQ : public BaseTask {
 public:
  static constexpr ppc::task::TypeOfTask GetStaticTypeOfTask() {
    return ppc::task::TypeOfTask::kSEQ;
  }
  explicit KamaletdinovRBitwiseIntSEQ(const InType &in);

 private:
  bool ValidationImpl() override;
  bool PreProcessingImpl() override;
  bool RunImpl() override;
  bool PostProcessingImpl() override;

  std::vector<int> data_;
};

}  // namespace kamaletdinov_r_bitwise_int_seq
