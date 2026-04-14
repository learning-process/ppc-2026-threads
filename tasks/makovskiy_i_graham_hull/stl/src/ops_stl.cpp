#include "makovskiy_i_graham_hull/stl/include/ops_stl.hpp"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <future>
#include <thread>
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

bool IsBetterMin(const Point &candidate, const Point &current_min) {
  if (candidate.y < current_min.y - 1e-9) {
    return true;
  }
  return (std::abs(candidate.y - current_min.y) <= 1e-9) && (candidate.x < current_min.x);
}

size_t FindMinPointIndexSTL(const std::vector<Point> &points) {
  size_t n = points.size();
  unsigned int num_threads = std::thread::hardware_concurrency();
  if (num_threads == 0) {
    num_threads = 4;
  }
  if (n < 1000) {
    num_threads = 1;
  }

  std::vector<std::future<size_t>> futures;
  size_t chunk = (n + num_threads - 1) / num_threads;

  auto worker = [&points](size_t start, size_t end) {
    size_t local_min = start;
    for (size_t j = start + 1; j < end; ++j) {
      if (IsBetterMin(points[j], points[local_min])) {
        local_min = j;
      }
    }
    return local_min;
  };

  for (unsigned int i = 0; i < num_threads; ++i) {
    size_t start = i * chunk;
    size_t end = std::min(start + chunk, n);
    if (start >= n) {
      break;
    }
    futures.push_back(std::async(std::launch::async, worker, start, end));
  }

  size_t min_idx = futures[0].get();
  for (size_t i = 1; i < futures.size(); ++i) {
    size_t local_min = futures[i].get();
    if (IsBetterMin(points[local_min], points[min_idx])) {
      min_idx = local_min;
    }
  }

  return min_idx;
}

template <typename RandomIt, typename Compare>
void StlParallelSortSub(RandomIt first, RandomIt last, Compare comp) {
  if (last - first < 2048) {
    std::sort(first, last, comp);
    return;
  }
  auto pivot = *(first + ((last - first) / 2));
  RandomIt middle1 = std::partition(first, last, [pivot, comp](const auto &a) { return comp(a, pivot); });
  RandomIt middle2 = std::partition(middle1, last, [pivot, comp](const auto &a) { return !comp(pivot, a); });

  auto future1 = std::async(std::launch::async, [first, middle1, comp]() { std::sort(first, middle1, comp); });
  std::sort(middle2, last, comp);
  future1.wait();
}

template <typename RandomIt, typename Compare>
void StlParallelSort(RandomIt first, RandomIt last, Compare comp) {
  if (last - first < 2048) {
    std::sort(first, last, comp);
    return;
  }
  auto pivot = *(first + ((last - first) / 2));
  RandomIt middle1 = std::partition(first, last, [pivot, comp](const auto &a) { return comp(a, pivot); });
  RandomIt middle2 = std::partition(middle1, last, [pivot, comp](const auto &a) { return !comp(pivot, a); });

  auto future1 = std::async(std::launch::async, [first, middle1, comp]() { StlParallelSortSub(first, middle1, comp); });
  StlParallelSortSub(middle2, last, comp);
  future1.wait();
}

std::vector<Point> FilterPointsSTL(const std::vector<Point> &points, const Point &p0) {
  size_t n = points.size();

  if (n <= 2) {
    return points;
  }

  std::vector<uint8_t> keep(n, 1);
  unsigned int num_threads = std::thread::hardware_concurrency();

  if (num_threads == 0) {
    num_threads = 4;
  }
  if (n < 1000) {
    num_threads = 1;
  }

  size_t num_elements = n - 2;
  size_t chunk = (num_elements + num_threads - 1) / num_threads;
  std::vector<std::future<void>> futures;

  auto worker = [&points, &keep, &p0](size_t start, size_t end) {
    for (size_t j = start; j < end; ++j) {
      if (std::abs(CrossProduct(p0, points[j], points[j + 1])) < 1e-9) {
        keep[j] = 0;
      }
    }
  };

  for (unsigned int i = 0; i < num_threads; ++i) {
    size_t start = 1 + (i * chunk);
    size_t end = std::min(start + chunk, n - 1);

    if (start >= n - 1) {
      break;
    }

    futures.push_back(std::async(std::launch::async, worker, start, end));
  }

  for (auto &f : futures) {
    f.wait();
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

ConvexHullGrahamSTL::ConvexHullGrahamSTL(const InType &in) {
  SetTypeOfTask(GetStaticTypeOfTask());
  GetInput() = in;
}

bool ConvexHullGrahamSTL::ValidationImpl() {
  return true;
}

bool ConvexHullGrahamSTL::PreProcessingImpl() {
  return true;
}

bool ConvexHullGrahamSTL::RunImpl() {
  InType points = GetInput();
  if (points.size() < 3) {
    GetOutput() = points;
    return true;
  }

  size_t min_idx = FindMinPointIndexSTL(points);

  std::swap(points[0], points[min_idx]);
  Point p0 = points[0];

  auto comp = [p0](const Point &a, const Point &b) {
    double cp = CrossProduct(p0, a, b);
    if (std::abs(cp) < 1e-9) {
      return DistSq(p0, a) < DistSq(p0, b);
    }
    return cp > 0;
  };

  StlParallelSort(points.begin() + 1, points.end(), comp);

  InType filtered = FilterPointsSTL(points, p0);

  if (filtered.size() < 3) {
    GetOutput() = filtered;
    return true;
  }

  GetOutput() = BuildHull(filtered);
  return true;
}

bool ConvexHullGrahamSTL::PostProcessingImpl() {
  return true;
}

}  // namespace makovskiy_i_graham_hull
