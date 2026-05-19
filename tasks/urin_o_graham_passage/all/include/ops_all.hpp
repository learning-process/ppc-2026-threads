#pragma once

#include <vector>

#include "task/include/task.hpp"
#include "urin_o_graham_passage/common/include/common.hpp"

namespace urin_o_graham_passage {

class UrinOGrahamPassageALL : public BaseTask {
 public:
  static constexpr ppc::task::TypeOfTask GetStaticTypeOfTask() {
    return ppc::task::TypeOfTask::kALL;
  }
  explicit UrinOGrahamPassageALL(const InType &in);

 protected:
  bool ValidationImpl() override;
  bool PreProcessingImpl() override;
  bool RunImpl() override;
  bool PostProcessingImpl() override;

 private:
  static Point FindLowestPoint(const InType &points);
  static double PolarAngle(const Point &base, const Point &p);
  static int Orientation(const Point &p, const Point &q, const Point &r);
  static double DistanceSquared(const Point &p1, const Point &p2);

  // Параллельные версии (используют все доступные технологии)
  static Point FindLowestPointParallel(const InType &points);
  static std::vector<Point> PrepareOtherPointsParallel(const InType &points, const Point &p0);
  static bool AreAllCollinear(const Point &p0, const std::vector<Point> &points);
  static std::vector<Point> BuildConvexHull(const Point &p0, const std::vector<Point> &points);
};

}  // namespace urin_o_graham_passage
