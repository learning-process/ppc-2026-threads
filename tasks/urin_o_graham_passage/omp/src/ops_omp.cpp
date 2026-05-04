#include "urin_o_graham_passage/omp/include/ops_omp.hpp"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <vector>

#ifdef _OPENMP
#  include <omp.h>
#endif

#include "urin_o_graham_passage/common/include/common.hpp"

namespace urin_o_graham_passage {

UrinOGrahamPassageOMP::UrinOGrahamPassageOMP(const InType &in) {
  SetTypeOfTask(GetStaticTypeOfTask());
  GetInput() = in;
  GetOutput() = OutType();
}

bool UrinOGrahamPassageOMP::ValidationImpl() {
  const auto &points = GetInput();

  if (points.size() < 3) {
    return false;
  }

  const Point &first = points[0];
  for (size_t i = 1; i < points.size(); ++i) {
    if (points[i] != first) {
      return true;
    }
  }

  return false;
}

bool UrinOGrahamPassageOMP::PreProcessingImpl() {
  GetOutput().clear();
  return true;
}

Point UrinOGrahamPassageOMP::FindLowestPoint(const InType &points) {
  Point lowest = points[0];

  for (size_t i = 1; i < points.size(); ++i) {
    if (points[i].y < lowest.y - 1e-10 ||
        (std::abs(points[i].y - lowest.y) < 1e-10 && points[i].x < lowest.x - 1e-10)) {
      lowest = points[i];
    }
  }

  return lowest;
}

Point UrinOGrahamPassageOMP::FindLowestPointParallel(const InType &points) {
  int lowest_index = 0;
  int n = static_cast<int>(points.size());  // Сохраняем размер в int

#pragma omp parallel for default(none) shared(points, lowest_index, n)
  for (int i = 1; i < n; ++i) {  // Простое сравнение
#pragma omp critical
    {
      if (points[i].y < points[lowest_index].y - 1e-10 ||
          (std::abs(points[i].y - points[lowest_index].y) < 1e-10 && points[i].x < points[lowest_index].x - 1e-10)) {
        lowest_index = i;
      }
    }
  }

  return points[lowest_index];
}

double UrinOGrahamPassageOMP::PolarAngle(const Point &base, const Point &p) {
  double dx = p.x - base.x;
  double dy = p.y - base.y;

  if (std::abs(dx) < 1e-10 && std::abs(dy) < 1e-10) {
    return -1e10;
  }

  return std::atan2(dy, dx);
}

int UrinOGrahamPassageOMP::Orientation(const Point &p, const Point &q, const Point &r) {
  double val = ((q.x - p.x) * (r.y - p.y)) - ((q.y - p.y) * (r.x - p.x));

  if (std::abs(val) < 1e-10) {
    return 0;
  }
  return (val > 0) ? 1 : -1;
}

double UrinOGrahamPassageOMP::DistanceSquared(const Point &p1, const Point &p2) {
  double dx = p2.x - p1.x;
  double dy = p2.y - p1.y;
  return (dx * dx) + (dy * dy);
}

std::vector<Point> UrinOGrahamPassageOMP::PrepareOtherPoints(const InType &points, const Point &p0) {
  std::vector<Point> other_points;
  other_points.reserve(points.size() - 1);

  for (const auto &point : points) {
    if (point != p0) {
      other_points.push_back(point);
    }
  }

  std::ranges::sort(other_points, [&p0](const Point &a, const Point &b) {  // Исправлено
    double angle_a = PolarAngle(p0, a);
    double angle_b = PolarAngle(p0, b);

    if (std::abs(angle_a - angle_b) < 1e-10) {
      return DistanceSquared(p0, a) < DistanceSquared(p0, b);
    }
    return angle_a < angle_b;
  });

  return other_points;
}

std::vector<Point> UrinOGrahamPassageOMP::PrepareOtherPointsParallel(const InType &points, const Point &p0) {
  std::vector<Point> other_points;
  other_points.reserve(points.size() - 1);

#ifdef _OPENMP
  int n = static_cast<int>(points.size());  // Сохраняем размер в int

#  pragma omp parallel default(none) shared(points, p0, other_points, n)
  {
    std::vector<Point> local_points;
    local_points.reserve((n / omp_get_num_threads()) + 1);

#  pragma omp for nowait
    for (int i = 0; i < n; ++i) {  // Простое сравнение
      if (points[i] != p0) {
        local_points.push_back(points[i]);
      }
    }

#  pragma omp critical
    {
      other_points.insert(other_points.end(), local_points.begin(), local_points.end());
    }
  }
#else
  for (const auto &point : points) {
    if (point != p0) {
      other_points.push_back(point);
    }
  }
#endif

  std::ranges::sort(other_points, [&p0](const Point &a, const Point &b) {  // Исправлено
    double angle_a = PolarAngle(p0, a);
    double angle_b = PolarAngle(p0, b);

    if (std::abs(angle_a - angle_b) < 1e-10) {
      return DistanceSquared(p0, a) < DistanceSquared(p0, b);
    }
    return angle_a < angle_b;
  });

  return other_points;
}

bool UrinOGrahamPassageOMP::AreAllCollinear(const Point &p0, const std::vector<Point> &points) {
  bool all_collinear = true;
  int n = static_cast<int>(points.size());  // Сохраняем размер в int

#pragma omp parallel for default(none) shared(points, p0, all_collinear, n)
  for (int i = 1; i < n; ++i) {  // Простое сравнение
    if (Orientation(p0, points[0], points[i]) != 0) {
#pragma omp atomic write
      all_collinear = false;
    }
  }

  return all_collinear;
}

std::vector<Point> UrinOGrahamPassageOMP::BuildConvexHull(const Point &p0, const std::vector<Point> &points) {
  std::vector<Point> hull;
  hull.reserve(points.size() + 1);
  hull.push_back(p0);
  hull.push_back(points[0]);

  for (size_t i = 1; i < points.size(); ++i) {
    while (hull.size() >= 2) {
      const Point &p = hull[hull.size() - 2];
      const Point &q = hull.back();
      if (Orientation(p, q, points[i]) <= 0) {
        hull.pop_back();
      } else {
        break;
      }
    }
    hull.push_back(points[i]);
  }

  return hull;
}

bool UrinOGrahamPassageOMP::RunImpl() {
  const InType &points = GetInput();

  if (points.size() < 3) {
    return false;
  }

  Point p0 = FindLowestPointParallel(points);

  std::vector<Point> other_points = PrepareOtherPointsParallel(points, p0);
  if (other_points.empty()) {
    return false;
  }

  if (AreAllCollinear(p0, other_points)) {
    GetOutput() = {p0, other_points.back()};
    return true;
  }

  GetOutput() = BuildConvexHull(p0, other_points);
  return true;
}

bool UrinOGrahamPassageOMP::PostProcessingImpl() {
  return !GetOutput().empty();
}

}  // namespace urin_o_graham_passage
