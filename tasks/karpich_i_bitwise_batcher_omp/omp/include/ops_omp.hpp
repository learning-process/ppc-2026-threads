#pragma once

#include <vector>

#include "karpich_i_bitwise_batcher_omp/common/include/common.hpp"
#include "task/include/task.hpp"

namespace karpich_i_bitwise_batcher_omp {

class KarpichIBitwiseBatcherOMP : public BaseTask {
 public:
  static constexpr ppc::task::TypeOfTask GetStaticTypeOfTask() {
    return ppc::task::TypeOfTask::kOMP;
  }
  explicit KarpichIBitwiseBatcherOMP(const InType &in);

 private:
  bool ValidationImpl() override;
  bool PreProcessingImpl() override;
  bool RunImpl() override;
  bool PostProcessingImpl() override;

  std::vector<int> data_;
};

}  // namespace karpich_i_bitwise_batcher_omp
