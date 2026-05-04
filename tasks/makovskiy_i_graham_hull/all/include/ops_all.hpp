#pragma once

#include "makovskiy_i_graham_hull/common/include/common.hpp"
#include "task/include/task.hpp"

namespace makovskiy_i_graham_hull {

class ConvexHullGrahamALL : public BaseTask {
 public:
  static constexpr ppc::task::TypeOfTask GetStaticTypeOfTask() {
    return ppc::task::TypeOfTask::kALL;
  }

  explicit ConvexHullGrahamALL(const InType &in);

 protected:
  bool ValidationImpl() override;
  bool PreProcessingImpl() override;
  bool RunImpl() override;
  bool PostProcessingImpl() override;
};

}  // namespace makovskiy_i_graham_hull
