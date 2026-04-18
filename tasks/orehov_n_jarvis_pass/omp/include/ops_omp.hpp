#pragma once

#include <vector>

#include "orehov_n_jarvis_pass/common/include/common.hpp"
#include "task/include/task.hpp"

namespace orehov_n_jarvis_pass {

class OrehovNJarvisPassOMP : public ppc::task::Task<std::vector<Point>, std::vector<Point>> {
 public:
  static constexpr ppc::task::TypeOfTask GetStaticTypeOfTask() {
    return ppc::task::TypeOfTask::kOMP;
  }
  explicit OrehovNJarvisPassOMP(const std::vector<Point> &in);

 private:
  static double CheckLeft(Point a, Point b, Point c);
  static double DistanceSquared(Point a, Point b);
  [[nodiscard]] static Point FindFirstElem(const std::vector<Point> &input);
  [[nodiscard]] static Point FindNext(Point current, const std::vector<Point> &input);
  static void UpdateBestCandidate(Point current, const Point &candidate, Point &best, double orient);

  bool ValidationImpl() override;
  bool PreProcessingImpl() override;
  bool RunImpl() override;
  bool PostProcessingImpl() override;
};

}  // namespace orehov_n_jarvis_pass
