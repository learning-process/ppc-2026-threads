#pragma once

#include <vector>

#include "egorova_l_binary_convex_hull/common/include/common.hpp"
#include "task/include/task.hpp"

namespace egorova_l_binary_convex_hull {

class BinaryConvexHullSEQ : public BaseTask {
 public:
  static constexpr ppc::task::TypeOfTask GetStaticTypeOfTask() {
    return ppc::task::TypeOfTask::kSEQ;
  }
  explicit BinaryConvexHullSEQ(const InType &in);

  bool ValidationImpl() override;
  bool PreProcessingImpl() override;
  bool RunImpl() override;
  bool PostProcessingImpl() override;
};

}  // namespace egorova_l_binary_convex_hull
