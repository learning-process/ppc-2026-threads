#pragma once

#include "nalitov_d_dijkstras_algorithm_seq/common/include/common.hpp"
#include "task/include/task.hpp"

namespace nalitov_d_dijkstras_algorithm_seq {

class NalitovDDijkstrasAlgorithmSeq : public BaseTask {
 public:
  static constexpr ppc::task::TypeOfTask GetStaticTypeOfTask() {
    return ppc::task::TypeOfTask::kSEQ;
  }
  explicit NalitovDDijkstrasAlgorithmSeq(const InType &in);

 private:
  bool ValidationImpl() override;
  bool PreProcessingImpl() override;
  bool RunImpl() override;
  bool PostProcessingImpl() override;
};

}  // namespace nalitov_d_dijkstras_algorithm_seq
