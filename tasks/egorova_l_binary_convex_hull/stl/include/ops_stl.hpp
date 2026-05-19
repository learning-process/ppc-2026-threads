#pragma once

#include <vector>

#include "egorova_l_binary_convex_hull/common/include/common.hpp"
#include "task/include/task.hpp"

namespace egorova_l_binary_convex_hull {

class BinaryConvexHullSTL : public BaseTask {
 public:
  static constexpr ppc::task::TypeOfTask GetStaticTypeOfTask() {
    return ppc::task::TypeOfTask::kSTL;
  }
  explicit BinaryConvexHullSTL(const InType &in);

 private:
  bool ValidationImpl() override;
  bool PreProcessingImpl() override;
  bool RunImpl() override;
  bool PostProcessingImpl() override;

  std::vector<std::vector<Point>> components_;

  static std::vector<std::vector<Point>> FindComponents(const ImageData &img);

  static std::vector<Point> ConvexHull(std::vector<Point> points);
};

}  // namespace egorova_l_binary_convex_hull
