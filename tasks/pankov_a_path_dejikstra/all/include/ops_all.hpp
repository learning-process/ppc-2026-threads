#pragma once

#include "pankov_a_path_dejikstra/common/include/common.hpp"
#include "task/include/task.hpp"

namespace pankov_a_path_dejikstra {

class PankovAPathDejikstraALL : public BaseTask {
 public:
  static constexpr ppc::task::TypeOfTask GetStaticTypeOfTask() {
    return ppc::task::TypeOfTask::kALL;
  }
  explicit PankovAPathDejikstraALL(const InType &in);

 private:
  bool ValidationImpl() override;
  bool PreProcessingImpl() override;
  bool RunImpl() override;
  bool PostProcessingImpl() override;
};

}  // namespace pankov_a_path_dejikstra
