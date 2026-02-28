#include "urin_o_graham_passage/seq/include/ops_seq.hpp"

#include <algorithm>
#include <cmath>
#include <vector>

namespace urin_o_graham_passage {

UrinOGrahamPassageSEQ::UrinOGrahamPassageSEQ(const InType &in) {
  SetTypeOfTask(GetStaticTypeOfTask());
  GetInput() = in;
  GetOutput() = OutType();
}

bool UrinOGrahamPassageSEQ::ValidationImpl() {
  const auto &points = GetInput();

  if (points.size() < 3) {
    return false;
  }

  // Проверяем, что не все точки одинаковые
  const Point &first = points[0];
  bool all_same = true;
  for (size_t i = 1; i < points.size(); i++) {
    if (!(points[i] == first)) {
      all_same = false;
      break;
    }
  }

  if (all_same) {
    return false;
  }

  return true;
}

bool UrinOGrahamPassageSEQ::PreProcessingImpl() {
  GetOutput().clear();
  return true;
}

Point UrinOGrahamPassageSEQ::FindLowestPoint(const InType &points) {
  Point lowest = points[0];

  for (size_t i = 1; i < points.size(); i++) {
    if (points[i].y < lowest.y - 1e-10 ||
        (std::abs(points[i].y - lowest.y) < 1e-10 && points[i].x < lowest.x - 1e-10)) {
      lowest = points[i];
    }
  }

  return lowest;
}

double UrinOGrahamPassageSEQ::PolarAngle(const Point &base, const Point &p) {
  double dx = p.x - base.x;
  double dy = p.y - base.y;

  if (std::abs(dx) < 1e-10 && std::abs(dy) < 1e-10) {
    return -1e10;
  }

  return std::atan2(dy, dx);
}

int UrinOGrahamPassageSEQ::Orientation(const Point &p, const Point &q, const Point &r) {
  double val = (q.x - p.x) * (r.y - p.y) - (q.y - p.y) * (r.x - p.x);

  if (std::abs(val) < 1e-10) {
    return 0;
  }
  return (val > 0) ? 1 : -1;
}

double UrinOGrahamPassageSEQ::DistanceSquared(const Point &p1, const Point &p2) {
  double dx = p2.x - p1.x;
  double dy = p2.y - p1.y;
  return dx * dx + dy * dy;
}

bool UrinOGrahamPassageSEQ::RunImpl() {
  const InType &points = GetInput();

  if (points.size() < 3) {
    return false;
  }

  Point p0 = FindLowestPoint(points);

  std::vector<Point> other_points;
  for (const auto &point : points) {
    if (!(point == p0)) {
      other_points.push_back(point);
    }
  }

  // ИСПРАВЛЕНО: убрали захват this
  std::sort(other_points.begin(), other_points.end(), [&p0](const Point &a, const Point &b) {
    double dx1 = a.x - p0.x;
    double dy1 = a.y - p0.y;
    double dx2 = b.x - p0.x;
    double dy2 = b.y - p0.y;

    double angle_a = std::atan2(dy1, dx1);
    double angle_b = std::atan2(dy2, dx2);

    if (std::abs(angle_a - angle_b) < 1e-10) {
      // Если углы равны, ближе та, у которой меньше расстояние
      return (dx1 * dx1 + dy1 * dy1) < (dx2 * dx2 + dy2 * dy2);
    }
    return angle_a < angle_b;
  });

  // Проверяем, не все ли точки коллинеарны
  bool all_collinear = true;
  for (size_t i = 1; i < other_points.size(); i++) {
    if (Orientation(p0, other_points[0], other_points[i]) != 0) {
      all_collinear = false;
      break;
    }
  }

  if (all_collinear) {
    GetOutput() = {p0, other_points.back()};
    return true;
  }

  // Строим выпуклую оболочку
  std::vector<Point> hull;
  hull.push_back(p0);
  hull.push_back(other_points[0]);

  for (size_t i = 1; i < other_points.size(); i++) {
    while (hull.size() >= 2) {
      Point &p = hull[hull.size() - 2];
      Point &q = hull[hull.size() - 1];
      int orient = Orientation(p, q, other_points[i]);

      if (orient <= 0) {
        hull.pop_back();
      } else {
        break;
      }
    }
    hull.push_back(other_points[i]);
  }

  GetOutput() = hull;
  return true;
}

bool UrinOGrahamPassageSEQ::PostProcessingImpl() {
  return !GetOutput().empty();
}

}  // namespace urin_o_graham_passage
