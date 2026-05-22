#pragma once

#include <cstdint>
#include <vector>

#include "../../common/include/common.hpp"
#include "task/include/task.hpp"

namespace kaur_a_dijkstra_alg {

class KaurADijkstraAlgALL : public BaseTask {
 public:
  static constexpr ppc::task::TypeOfTask GetStaticTypeOfTask() {
    return ppc::task::TypeOfTask::kALL;
  }
  explicit KaurADijkstraAlgALL(const InType &in);

 private:
  bool ValidationImpl() override;
  bool PreProcessingImpl() override;
  bool RunImpl() override;
  bool PostProcessingImpl() override;

  std::vector<InType> dist_;
  std::vector<char> visited_;

  int rank_{};
  int size_{};
  int local_n_{};
  int start_v_{};
  int64_t total_sum_{};
};

}  // namespace kaur_a_dijkstra_alg
