#pragma once

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

  std::vector<InType> distances_;
  std::vector<char> processed_;
};

}  // namespace nalitov_d_dijkstras_algorithm
