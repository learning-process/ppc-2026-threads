#pragma once

#include <vector>

#include "../../common/include/common.hpp"
#include "task/include/task.hpp"

namespace kaur_a_dijkstra_alg {

class KaurADijkstraAlgSEQ : public BaseTask {
 public:
  static constexpr ppc::task::TypeOfTask GetStaticTypeOfTask() {
    return ppc::task::TypeOfTask::kSEQ;
  }
  explicit KaurADijkstraAlgSEQ(const InType &in);

 private:
  bool ValidationImpl() override;
  bool PreProcessingImpl() override;
  bool RunImpl() override;
  bool PostProcessingImpl() override;
  static InType FindMinDist(const std::vector<InType> &dist, const std::vector<bool> &visited);
  static void RelaxEdges(InType u, std::vector<InType> &dist, const std::vector<bool> &visited);
};

}  // namespace kaur_a_dijkstra_alg
