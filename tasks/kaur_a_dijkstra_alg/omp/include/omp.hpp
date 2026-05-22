#pragma once

#include <vector>

#include "../../common/include/common.hpp"
#include "task/include/task.hpp"

namespace kaur_a_dijkstra_alg {

class KaurADijkstraAlgOMP : public BaseTask {
 public:
  static constexpr ppc::task::TypeOfTask GetStaticTypeOfTask() {
    return ppc::task::TypeOfTask::kOMP;
  }
  explicit KaurADijkstraAlgOMP(const InType &in);

 private:
  bool ValidationImpl() override;
  bool PreProcessingImpl() override;
  bool RunImpl() override;
  bool PostProcessingImpl() override;

  std::vector<InType> dist_;
  std::vector<char> visited_;
};

}  // namespace kaur_a_dijkstra_alg
