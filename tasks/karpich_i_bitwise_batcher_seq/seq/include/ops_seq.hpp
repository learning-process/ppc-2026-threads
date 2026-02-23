#pragma once

#include <vector>

#include "karpich_i_bitwise_batcher_seq/common/include/common.hpp"
#include "task/include/task.hpp"

namespace karpich_i_bitwise_batcher_seq {

class KarpichIBitwiseBatcherSEQ : public BaseTask {
 public:
  static constexpr ppc::task::TypeOfTask GetStaticTypeOfTask() {
    return ppc::task::TypeOfTask::kSEQ;
  }
  explicit KarpichIBitwiseBatcherSEQ(const InType &in);

 private:
  bool ValidationImpl() override;
  bool PreProcessingImpl() override;
  bool RunImpl() override;
  bool PostProcessingImpl() override;

  std::vector<int> data_;
};

}  // namespace karpich_i_bitwise_batcher_seq
