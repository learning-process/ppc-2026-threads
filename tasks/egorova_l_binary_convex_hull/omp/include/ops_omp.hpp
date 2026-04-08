#pragma once

#include <cstdint>
#include <vector>

#include "egorova_l_binary_convex_hull/common/include/common.hpp"
#include "task/include/task.hpp"

namespace egorova_l_binary_convex_hull {

class BinaryConvexHullOMP : public BaseTask {
 public:
  static constexpr ppc::task::TypeOfTask GetStaticTypeOfTask() {
    return ppc::task::TypeOfTask::kOMP;
  }
  explicit BinaryConvexHullOMP(const InType &in);

  bool ValidationImpl() override;
  bool PreProcessingImpl() override;
  bool RunImpl() override;
  bool PostProcessingImpl() override;

 private:
  // Вспомогательные функции
  static std::vector<std::vector<Point>> FindComponents(const std::vector<uint8_t> &image, int width, int height);
  static void ProcessComponent(const std::vector<uint8_t> &image, int width, int height, int start_x, int start_y,
                               int label, std::vector<int> &labels, std::vector<Point> &component_points);
  static void BuildConvexHull(std::vector<Point> &points, std::vector<Point> &hull);
  static int64_t CrossProduct(const Point &a, const Point &b, const Point &c);
};

}  // namespace egorova_l_binary_convex_hull
