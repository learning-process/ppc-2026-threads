#include "makovskiy_i_graham_hull/seq/include/ops_seq.hpp"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <vector>

#include "makovskiy_i_graham_hull/common/include/common.hpp"

namespace makovskiy_i_graham_hull {

namespace {

double CrossProduct(const Point &o, const Point &a, const Point &b) {
  return ((a.x - o.x) * (b.y - o.y)) - ((a.y - o.y) * (b.x - o.x));
}

double DistSq(const Point &a, const Point &b) {
  return ((a.x - b.x) * (a.x - b.x)) + ((a.y - b.y) * (a.y - b.y));
}

size_t FindMinPointIndex(const std::vector<Point> &points) {
  size_t min_idx = 0;
  for (size_t i = 1; i < points.size(); ++i) {
    if (points[i].y < points[min_idx].y - 1e-9 ||
        (std::abs(points[i].y - points[min_idx].y) <= 1e-9 && points[i].x < points[min_idx].x)) {
      min_idx = i;
    }
  }
  return min_idx;
}

std::vector<Point> FilterPoints(const std::vector<Point> &points, const Point &p0) {
  std::vector<Point> filtered;
  filtered.push_back(points[0]);
  for (size_t i = 1; i < points.size(); ++i) {
    while (i < points.size() - 1 && std::abs(CrossProduct(p0, points[i], points[i + 1])) < 1e-9) {
      i++;
    }
    filtered.push_back(points[i]);
  }
  return filtered;
}

std::vector<Point> BuildHull(const std::vector<Point> &filtered) {
  std::vector<Point> hull;
  hull.push_back(filtered[0]);
  hull.push_back(filtered[1]);
  hull.push_back(filtered[2]);

  for (size_t i = 3; i < filtered.size(); ++i) {
    while (hull.size() > 1 && CrossProduct(hull[hull.size() - 2], hull.back(), filtered[i]) <= 1e-9) {
      hull.pop_back();
    }
    hull.push_back(filtered[i]);
  }
  return hull;
}

}  // namespace

ConvexHullGrahamSEQ::ConvexHullGrahamSEQ(const InType &in) {
  SetTypeOfTask(GetStaticTypeOfTask());
  GetInput() = in;
}

bool ConvexHullGrahamSEQ::ValidationImpl() {
  return true;
}

bool ConvexHullGrahamSEQ::PreProcessingImpl() {
  return true;
}

bool ConvexHullGrahamSEQ::RunImpl() {
  InType points = GetInput();
  if (points.size() < 3) {
    GetOutput() = points;
    return true;
  }

  size_t min_idx = FindMinPointIndex(points);

  std::swap(points[0], points[min_idx]);
  Point p0 = points[0];

  std::sort(points.begin() + 1, points.end(), [&p0](const Point &a, const Point &b) {
    double cp = CrossProduct(p0, a, b);
    if (std::abs(cp) < 1e-9) {
      return DistSq(p0, a) < DistSq(p0, b);
    }
    return cp > 0;
  });

  InType filtered = FilterPoints(points, p0);

  if (filtered.size() < 3) {
    GetOutput() = filtered;
    return true;
  }

  GetOutput() = BuildHull(filtered);
  return true;
}

bool ConvexHullGrahamSEQ::PostProcessingImpl() {
  return true;
}

}  // namespace makovskiy_i_graham_hull
