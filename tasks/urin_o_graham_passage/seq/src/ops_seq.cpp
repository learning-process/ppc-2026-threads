#include "urin_o_graham_passage/seq/include/ops_seq.hpp"

#include <algorithm>
#include <cmath>
#include <cstddef>  // Добавлено для size_t
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

  const Point &first = points[0];
  for (size_t i = 1; i < points.size(); ++i) {
    if (points[i] != first) {
      return true;  // Нашли разную точку - валидация пройдена
    }
  }

  return false;  // Все точки одинаковые
}

bool UrinOGrahamPassageSEQ::PreProcessingImpl() {
  GetOutput().clear();
  return true;
}

Point UrinOGrahamPassageSEQ::FindLowestPoint(const InType &points) {
  Point lowest = points[0];

  for (size_t i = 1; i < points.size(); ++i) {
    if (points[i].y < lowest.y - 1e-10 ||
        (std::abs(points[i].y - lowest.y) < 1e-10 && points[i].x < lowest.x - 1e-10)) {
      lowest = points[i];
    }
  }

  return lowest;
}

double UrinOGrahamPassageSEQ::PolarAngle(const Point &base, const Point &p) {
  const double dx = p.x - base.x;
  const double dy = p.y - base.y;

  if (std::abs(dx) < 1e-10 && std::abs(dy) < 1e-10) {
    return -1e10;
  }

  return std::atan2(dy, dx);
}

int UrinOGrahamPassageSEQ::Orientation(const Point &p, const Point &q, const Point &r) {
  const double val = ((q.x - p.x) * (r.y - p.y)) - ((q.y - p.y) * (r.x - p.x));

  if (std::abs(val) < 1e-10) {
    return 0;
  }
  return (val > 0) ? 1 : -1;
}

double UrinOGrahamPassageSEQ::DistanceSquared(const Point &p1, const Point &p2) {
  const double dx = p2.x - p1.x;
  const double dy = p2.y - p1.y;
  return (dx * dx) + (dy * dy);  // Добавлены скобки для читаемости
}

bool UrinOGrahamPassageSEQ::RunImpl() {
  const InType &points = GetInput();

  if (points.size() < 3) {
    return false;
  }

  const Point p0 = FindLowestPoint(points);

  std::vector<Point> other_points;
  other_points.reserve(points.size() - 1);

  for (const auto &point : points) {
    if (point != p0) {
      other_points.push_back(point);
    }
  }

  std::sort(other_points.begin(), other_points.end(), [&p0](const Point &a, const Point &b) {
    const double dx1 = a.x - p0.x;
    const double dy1 = a.y - p0.y;
    const double dx2 = b.x - p0.x;
    const double dy2 = b.y - p0.y;

    const double angle_a = std::atan2(dy1, dx1);
    const double angle_b = std::atan2(dy2, dx2);

    if (std::abs(angle_a - angle_b) < 1e-10) {
      const double dist_a = (dx1 * dx1) + (dy1 * dy1);
      const double dist_b = (dx2 * dx2) + (dy2 * dy2);
      return dist_a < dist_b;
    }
    return angle_a < angle_b;
  });

  bool all_collinear = true;
  for (size_t i = 1; i < other_points.size(); ++i) {
    if (Orientation(p0, other_points[0], other_points[i]) != 0) {
      all_collinear = false;
      break;
    }
  }

  if (all_collinear) {
    GetOutput() = {p0, other_points.back()};
    return true;
  }

  std::vector<Point> hull;
  hull.reserve(other_points.size() + 1);
  hull.push_back(p0);
  hull.push_back(other_points[0]);

  for (size_t i = 1; i < other_points.size(); ++i) {
    while (hull.size() >= 2) {
      const Point &p = hull[hull.size() - 2];
      const Point &q = hull.back();
      const int orient = Orientation(p, q, other_points[i]);

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
