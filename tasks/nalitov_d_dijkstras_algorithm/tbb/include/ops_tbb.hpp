#pragma once

#include <utility>
#include <vector>

#include "nalitov_d_dijkstras_algorithm/common/include/common.hpp"
#include "task/include/task.hpp"

namespace nalitov_d_dijkstras_algorithm {

class NalitovDDijkstrasAlgorithmTBB : public BaseTask {
 public:
  static constexpr ppc::task::TypeOfTask GetStaticTypeOfTask() {
    return ppc::task::TypeOfTask::kTBB;
  }
  explicit NalitovDDijkstrasAlgorithmTBB(const InType &in);

 private:
  bool ValidationImpl() override;
  bool PreProcessingImpl() override;
  bool RunImpl() override;
  bool PostProcessingImpl() override;

  using OutgoingTable = std::vector<std::vector<std::pair<NodeId, Cost>>>;

  OutgoingTable graph_;
  std::vector<Cost> dist_;
  std::vector<char> visited_;
};

}  // namespace nalitov_d_dijkstras_algorithm
