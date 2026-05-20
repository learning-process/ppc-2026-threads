#pragma once

#include <utility>
#include <vector>

#include "nalitov_d_dijkstras_algorithm/common/include/common.hpp"
#include "task/include/task.hpp"

namespace nalitov_d_dijkstras_algorithm {

class NalitovDDijkstrasAlgorithmALL : public BaseTask {
 public:
  static constexpr ppc::task::TypeOfTask GetStaticTypeOfTask() {
    return ppc::task::TypeOfTask::kALL;
  }
  explicit NalitovDDijkstrasAlgorithmALL(const InType &in);

 private:
  bool ValidationImpl() override;
  bool PreProcessingImpl() override;
  bool RunImpl() override;
  bool PostProcessingImpl() override;

  using OutgoingTable = std::vector<std::vector<std::pair<NodeId, Cost>>>;
  OutgoingTable graph_;

  int n_{};
  int source_{};
  int rank_{};
  int size_{};
  int local_start_{};
  int local_count_{};
  std::vector<Cost> local_dist_;
  std::vector<char> local_visited_;
};

}  // namespace nalitov_d_dijkstras_algorithm
