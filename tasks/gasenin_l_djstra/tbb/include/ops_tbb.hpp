#pragma once

#include <limits>
#include <utility>
#include <vector>

#include "gasenin_l_djstra/common/include/common.hpp"
#include "task/include/task.hpp"

namespace gasenin_l_djstra {

class GaseninLDjstraTBB : public BaseTask {
 public:
  static constexpr ppc::task::TypeOfTask GetStaticTypeOfTask() {
    return ppc::task::TypeOfTask::kTBB;
  }
  explicit GaseninLDjstraTBB(const InType &in);

 private:
  bool ValidationImpl() override;
  bool PreProcessingImpl() override;
  bool RunImpl() override;
  bool PostProcessingImpl() override;

  int vertex_count_ = 0;
  std::vector<std::vector<std::pair<int, int>>> adj_;
  std::vector<int> dist_;
  static constexpr int kInf = std::numeric_limits<int>::max() / 2;
};

}  // namespace gasenin_l_djstra
