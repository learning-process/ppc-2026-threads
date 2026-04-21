#pragma once

#include "task/include/task.hpp"
#include "tsarkov_k_jarvis_convex_hull/common/include/common.hpp"

namespace tsarkov_k_jarvis_convex_hull {

class TsarkovKJarvisConvexHullTBB : public BaseTask {
 public:
  static constexpr ppc::task::TypeOfTask GetStaticTypeOfTask() {
    return ppc::task::TypeOfTask::kTBB;
  }

  explicit TsarkovKJarvisConvexHullTBB(const InType &input_points);

 private:
  bool ValidationImpl() override;
  bool PreProcessingImpl() override;
  bool RunImpl() override;
  bool PostProcessingImpl() override;
};

}  // namespace tsarkov_k_jarvis_convex_hull
