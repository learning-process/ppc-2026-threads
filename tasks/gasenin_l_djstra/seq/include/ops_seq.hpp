#pragma once

#include <vector>

#include "gasenin_l_djstra/common/include/common.hpp"
#include "task/include/task.hpp"

namespace gasenin_l_djstra {

class GaseninLDjstraSEQ : public BaseTask {
 public:
  static constexpr ppc::task::TypeOfTask GetStaticTypeOfTask() {
    return ppc::task::TypeOfTask::kSEQ;
  }
  explicit GaseninLDjstraSEQ(const InType &in);

 private:
  bool ValidationImpl() override;
  bool PreProcessingImpl() override;
  bool RunImpl() override;
  bool PostProcessingImpl() override;
  static InType FindMinDist(const std::vector<InType> &dist, const std::vector<bool> &visited);
  static void RelaxEdges(InType u, std::vector<InType> &dist, const std::vector<bool> &visited);
};

}  // namespace gasenin_l_djstra
