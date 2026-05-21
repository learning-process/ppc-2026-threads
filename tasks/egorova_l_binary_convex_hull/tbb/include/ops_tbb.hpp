#pragma once

#include "egorova_l_binary_convex_hull/common/include/common.hpp"
#include "task/include/task.hpp"

namespace egorova_l_binary_convex_hull {

class BinaryConvexHullTBB : public BaseTask {
 public:
  static constexpr ppc::task::TypeOfTask GetStaticTypeOfTask() {
    return ppc::task::TypeOfTask::kTBB;
  }
  explicit BinaryConvexHullTBB(const InType &in);

  bool ValidationImpl() override;
  bool PreProcessingImpl() override;
  bool RunImpl() override;
  bool PostProcessingImpl() override;
};

}  // namespace egorova_l_binary_convex_hull
