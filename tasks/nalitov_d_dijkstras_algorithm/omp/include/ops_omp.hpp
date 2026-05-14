#pragma once

#include <utility>
#include <vector>

#include "nalitov_d_dijkstras_algorithm/common/include/common.hpp"
#include "task/include/task.hpp"

namespace nalitov_d_dijkstras_algorithm {

class NalitovDDijkstrasAlgorithmOmp : public BaseTask {
 public:
  static constexpr ppc::task::TypeOfTask GetStaticTypeOfTask() {
    return ppc::task::TypeOfTask::kOMP;
  }
  explicit NalitovDDijkstrasAlgorithmOmp(const InType &in);

 private:
  bool ValidationImpl() override;
  bool PreProcessingImpl() override;
  bool RunImpl() override;
  bool PostProcessingImpl() override;

  using OutgoingTable = std::vector<std::vector<std::pair<NodeId, Cost>>>;
  OutgoingTable graph_;
};

}  // namespace nalitov_d_dijkstras_algorithm
