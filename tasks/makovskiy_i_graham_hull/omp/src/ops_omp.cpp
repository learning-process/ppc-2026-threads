#include "makovskiy_i_graham_hull/omp/include/ops_omp.hpp"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <cstdint>
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

size_t FindMinPointIndexOMP(const std::vector<Point> &points) {
  size_t min_idx = 0;
  size_t n = points.size();
#pragma omp parallel default(none) shared(points, n, min_idx)
  {
    size_t local_min = 0;
#pragma omp for
    for (size_t i = 1; i < n; ++i) {
      if (points[i].y < points[local_min].y - 1e-9 ||
          (std::abs(points[i].y - points[local_min].y) <= 1e-9 && points[i].x < points[local_min].x)) {
        local_min = i;
      }
    }
#pragma omp critical
    {
      if (points[local_min].y < points[min_idx].y - 1e-9 ||
          (std::abs(points[local_min].y - points[min_idx].y) <= 1e-9 && points[local_min].x < points[min_idx].x)) {
        min_idx = local_min;
      }
    }
  }
  return min_idx;
}

template <typename RandomIt, typename Compare>
void OmpQuickSort(RandomIt first, RandomIt last, Compare comp) {
  if (last - first < 2048) {
    std::sort(first, last, comp);
    return;
  }
  auto pivot = *(first + ((last - first) / 2));
  RandomIt middle1 = std::partition(first, last, [pivot, comp](const auto &a) { return comp(a, pivot); });
  RandomIt middle2 = std::partition(middle1, last, [pivot, comp](const auto &a) { return !comp(pivot, a); });

#pragma omp task default(none) firstprivate(first, middle1, comp)
  OmpQuickSort(first, middle1, comp);

#pragma omp task default(none) firstprivate(middle2, last, comp)
  OmpQuickSort(middle2, last, comp);

#pragma omp taskwait
}

std::vector<Point> FilterPointsOMP(const std::vector<Point> &points, const Point &p0) {
  size_t n = points.size();
  std::vector<uint8_t> keep(n, 1);

#pragma omp parallel for default(none) shared(n, points, keep, p0)
  for (size_t i = 1; i < n - 1; ++i) {
    if (std::abs(CrossProduct(p0, points[i], points[i + 1])) < 1e-9) {
      keep[i] = 0;
    }
  }

  std::vector<Point> filtered;
  filtered.reserve(n);
  for (size_t i = 0; i < n; ++i) {
    if (keep[i] != 0) {
      filtered.push_back(points[i]);
    }
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

ConvexHullGrahamOMP::ConvexHullGrahamOMP(const InType &in) {
  SetTypeOfTask(GetStaticTypeOfTask());
  GetInput() = in;
}

bool ConvexHullGrahamOMP::ValidationImpl() {
  return true;
}

bool ConvexHullGrahamOMP::PreProcessingImpl() {
  return true;
}

bool ConvexHullGrahamOMP::RunImpl() {
  InType points = GetInput();
  if (points.size() < 3) {
    GetOutput() = points;
    return true;
  }

  size_t min_idx = FindMinPointIndexOMP(points);

  std::swap(points[0], points[min_idx]);
  Point p0 = points[0];

  auto comp = [p0](const Point &a, const Point &b) {
    double cp = CrossProduct(p0, a, b);
    if (std::abs(cp) < 1e-9) {
      return DistSq(p0, a) < DistSq(p0, b);
    }
    return cp > 0;
  };

#pragma omp parallel default(none) shared(points, comp)
  {
#pragma omp single nowait
    OmpQuickSort(points.begin() + 1, points.end(), comp);
  }

  InType filtered = FilterPointsOMP(points, p0);

  if (filtered.size() < 3) {
    GetOutput() = filtered;
    return true;
  }

  GetOutput() = BuildHull(filtered);
  return true;
}

bool ConvexHullGrahamOMP::PostProcessingImpl() {
  return true;
}

}  // namespace makovskiy_i_graham_hull
