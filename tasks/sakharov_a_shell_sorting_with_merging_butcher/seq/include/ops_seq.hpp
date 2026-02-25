#pragma once

#include <vector>

#include "sakharov_a_shell_sorting_with_merging_butcher/common/include/common.hpp"
#include "task/include/task.hpp"

namespace sakharov_a_shell_sorting_with_merging_butcher {

class SakharovAShellButcherSEQ : public BaseTask {
 public:
  static constexpr ppc::task::TypeOfTask GetStaticTypeOfTask() {
    return ppc::task::TypeOfTask::kSEQ;
  }
  explicit SakharovAShellButcherSEQ(const InType &in);

 private:
  bool ValidationImpl() override;
  bool PreProcessingImpl() override;
  bool RunImpl() override;
  bool PostProcessingImpl() override;
  static InType FindMinDist(const std::vector<InType> &dist, const std::vector<bool> &visited);
  static void RelaxEdges(InType u, std::vector<InType> &dist, const std::vector<bool> &visited);
};

}  // namespace sakharov_a_shell_sorting_with_merging_butcher
