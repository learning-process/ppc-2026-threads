#pragma once

#include "kamalagin_a_binary_image_convex_hull_omp/common/include/common.hpp"
#include "task/include/task.hpp"

namespace kamalagin_a_binary_image_convex_hull_omp {

class KamalaginABinaryImageConvexHullOMP : public BaseTask {
 public:
  static constexpr ppc::task::TypeOfTask GetStaticTypeOfTask() {
    return ppc::task::TypeOfTask::kOMP;
  }
  explicit KamalaginABinaryImageConvexHullOMP(const InType &in);

 private:
  bool ValidationImpl() override;
  bool PreProcessingImpl() override;
  bool RunImpl() override;
  bool PostProcessingImpl() override;
};

}  // namespace kamalagin_a_binary_image_convex_hull_omp
