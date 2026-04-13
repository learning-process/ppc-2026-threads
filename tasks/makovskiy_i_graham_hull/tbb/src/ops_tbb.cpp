#include "makovskiy_i_graham_hull/tbb/include/ops_tbb.hpp"

#include <tbb/tbb.h>

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

size_t FindMinPointIndexTBB(const std::vector<Point> &points) {
  return tbb::parallel_reduce(tbb::blocked_range<size_t>(1, points.size()), static_cast<size_t>(0),
                              [&points](const tbb::blocked_range<size_t> &r, size_t local_min) {
    for (size_t i = r.begin(); i != r.end(); ++i) {
      if (points[i].y < points[local_min].y - 1e-9 ||
          (std::abs(points[i].y - points[local_min].y) <= 1e-9 && points[i].x < points[local_min].x)) {
        local_min = i;
      }
    }
    return local_min;
  }, [&points](size_t a, size_t b) {
    if (points[a].y < points[b].y - 1e-9 ||
        (std::abs(points[a].y - points[b].y) <= 1e-9 && points[a].x < points[b].x)) {
      return a;
    }
    return b;
  });
}

template <typename RandomIt, typename Compare>
void TbbQuickSort(RandomIt first, RandomIt last, Compare comp) {
  if (last - first < 2048) {
    std::sort(first, last, comp);
    return;
  }
  auto pivot = *(first + ((last - first) / 2));
  RandomIt middle1 = std::partition(first, last, [pivot, comp](const auto &a) { return comp(a, pivot); });
  RandomIt middle2 = std::partition(middle1, last, [pivot, comp](const auto &a) { return !comp(pivot, a); });

  tbb::task_group tg;
  tg.run([first, middle1, comp]() { TbbQuickSort(first, middle1, comp); });
  tg.run([middle2, last, comp]() { TbbQuickSort(middle2, last, comp); });
  tg.wait();
}

std::vector<Point> FilterPointsTBB(const std::vector<Point> &points, const Point &p0) {
  size_t n = points.size();
  std::vector<uint8_t> keep(n, 1);

  tbb::parallel_for(tbb::blocked_range<size_t>(1, n - 1), [&](const tbb::blocked_range<size_t> &r) {
    for (size_t i = r.begin(); i != r.end(); ++i) {
      if (std::abs(CrossProduct(p0, points[i], points[i + 1])) < 1e-9) {
        keep[i] = 0;
      }
    }
  });

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

ConvexHullGrahamTBB::ConvexHullGrahamTBB(const InType &in) {
  SetTypeOfTask(GetStaticTypeOfTask());
  GetInput() = in;
}

bool ConvexHullGrahamTBB::ValidationImpl() {
  return true;
}

bool ConvexHullGrahamTBB::PreProcessingImpl() {
  return true;
}

bool ConvexHullGrahamTBB::RunImpl() {
  InType points = GetInput();
  if (points.size() < 3) {
    GetOutput() = points;
    return true;
  }

  size_t min_idx = FindMinPointIndexTBB(points);

  std::swap(points[0], points[min_idx]);
  Point p0 = points[0];

  auto comp = [p0](const Point &a, const Point &b) {
    double cp = CrossProduct(p0, a, b);
    if (std::abs(cp) < 1e-9) {
      return DistSq(p0, a) < DistSq(p0, b);
    }
    return cp > 0;
  };

  TbbQuickSort(points.begin() + 1, points.end(), comp);

  InType filtered = FilterPointsTBB(points, p0);

  if (filtered.size() < 3) {
    GetOutput() = filtered;
    return true;
  }

  GetOutput() = BuildHull(filtered);
  return true;
}

bool ConvexHullGrahamTBB::PostProcessingImpl() {
  return true;
}

}  // namespace makovskiy_i_graham_hull
