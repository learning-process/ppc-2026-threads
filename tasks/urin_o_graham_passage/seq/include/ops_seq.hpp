#pragma once

#include <vector>

#include "task/include/task.hpp"
#include "urin_o_graham_passage/common/include/common.hpp"

namespace urin_o_graham_passage {

class UrinOGrahamPassageSEQ : public BaseTask {
 public:
  static constexpr ppc::task::TypeOfTask GetStaticTypeOfTask() {
    return ppc::task::TypeOfTask::kSEQ;
  }

  explicit UrinOGrahamPassageSEQ(const InType &in);

  [[nodiscard]] static Point FindLowestPoint(const InType &points);
  [[nodiscard]] static double PolarAngle(const Point &base, const Point &p);
  [[nodiscard]] static int Orientation(const Point &p, const Point &q, const Point &r);
  [[nodiscard]] static double DistanceSquared(const Point &p1, const Point &p2);

 protected:
  bool ValidationImpl() override;
  bool PreProcessingImpl() override;
  bool RunImpl() override;
  bool PostProcessingImpl() override;

 private:
  // Статические вспомогательные функции
  [[nodiscard]] static std::vector<Point> PrepareOtherPoints(const InType &points, const Point &p0);
  [[nodiscard]] static bool AreAllCollinear(const Point &p0, const std::vector<Point> &points);
  [[nodiscard]] static std::vector<Point> BuildConvexHull(const Point &p0, const std::vector<Point> &points);
};

}  // namespace urin_o_graham_passage
