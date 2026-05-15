#include "urin_o_graham_passage/stl/include/ops_stl.hpp"

#include <algorithm>
#include <atomic>
#include <cmath>
#include <cstddef>
#include <thread>
#include <vector>

#include "urin_o_graham_passage/common/include/common.hpp"

namespace urin_o_graham_passage {

UrinOGrahamPassageSTL::UrinOGrahamPassageSTL(const InType &in) {
  SetTypeOfTask(GetStaticTypeOfTask());
  GetInput() = in;
  GetOutput() = OutType();
}

bool UrinOGrahamPassageSTL::ValidationImpl() {
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

bool UrinOGrahamPassageSTL::PreProcessingImpl() {
  GetOutput().clear();
  return true;
}
//.
Point UrinOGrahamPassageSTL::FindLocalMinimum(const InType &points, size_t start, size_t end) {
  Point local_min = points[start];
  for (size_t i = start + 1; i < end; ++i) {
    if (points[i].y < local_min.y - 1e-10 ||
        (std::abs(points[i].y - local_min.y) < 1e-10 && points[i].x < local_min.x - 1e-10)) {
      local_min = points[i];
    }
  }
  return local_min;
}

Point UrinOGrahamPassageSTL::CombineMinimums(const std::vector<Point> &mins) {
  Point global_min = mins[0];
  for (size_t i = 1; i < mins.size(); ++i) {
    if (mins[i].y < global_min.y - 1e-10 ||
        (std::abs(mins[i].y - global_min.y) < 1e-10 && mins[i].x < global_min.x - 1e-10)) {
      global_min = mins[i];
    }
  }
  return global_min;
}

Point UrinOGrahamPassageSTL::FindLowestPoint(const InType &points) {
  Point lowest = points[0];

  for (size_t i = 1; i < points.size(); ++i) {
    if (points[i].y < lowest.y - 1e-10 ||
        (std::abs(points[i].y - lowest.y) < 1e-10 && points[i].x < lowest.x - 1e-10)) {
      lowest = points[i];
    }
  }

  return lowest;
}

Point UrinOGrahamPassageSTL::FindLowestPointParallel(const InType &points) {
  unsigned int num_threads = std::thread::hardware_concurrency();
  if (num_threads == 0) {
    num_threads = 2;
  }

  std::vector<std::thread> threads;
  std::vector<Point> local_mins(num_threads, points[0]);
  size_t block_size = points.size() / num_threads;

  for (unsigned int thread_idx = 0; thread_idx < num_threads; ++thread_idx) {
    size_t start = thread_idx * block_size;
    size_t end = (thread_idx == num_threads - 1) ? points.size() : start + block_size;

    threads.emplace_back([&points, start, end, &local_mins, thread_idx]() {
      local_mins[thread_idx] = FindLocalMinimum(points, start, end);
    });
  }

  for (auto &th : threads) {
    th.join();
  }

  return CombineMinimums(local_mins);
}

double UrinOGrahamPassageSTL::PolarAngle(const Point &base, const Point &p) {
  double dx = p.x - base.x;
  double dy = p.y - base.y;

  if (std::abs(dx) < 1e-10 && std::abs(dy) < 1e-10) {
    return -1e10;
  }

  return std::atan2(dy, dx);
}

int UrinOGrahamPassageSTL::Orientation(const Point &p, const Point &q, const Point &r) {
  double val = ((q.x - p.x) * (r.y - p.y)) - ((q.y - p.y) * (r.x - p.x));

  if (std::abs(val) < 1e-10) {
    return 0;
  }
  return (val > 0) ? 1 : -1;
}

double UrinOGrahamPassageSTL::DistanceSquared(const Point &p1, const Point &p2) {
  double dx = p2.x - p1.x;
  double dy = p2.y - p1.y;
  return (dx * dx) + (dy * dy);
}

std::vector<Point> UrinOGrahamPassageSTL::PrepareOtherPointsParallel(const InType &points, const Point &p0) {
  unsigned int num_threads = std::thread::hardware_concurrency();
  if (num_threads == 0) {
    num_threads = 2;
  }

  std::vector<std::thread> threads;
  std::vector<std::vector<Point>> local_points(num_threads);
  size_t block_size = points.size() / num_threads;

  for (unsigned int thread_idx = 0; thread_idx < num_threads; ++thread_idx) {
    size_t start = thread_idx * block_size;
    size_t end = (thread_idx == num_threads - 1) ? points.size() : start + block_size;

    threads.emplace_back([&points, &p0, start, end, &local_points, thread_idx]() {
      for (size_t i = start; i < end; ++i) {
        if (points[i] != p0) {
          local_points[thread_idx].push_back(points[i]);
        }
      }
    });
  }

  for (auto &th : threads) {
    th.join();
  }

  std::vector<Point> other_points;
  for (const auto &vec : local_points) {
    other_points.insert(other_points.end(), vec.begin(), vec.end());
  }

  std::ranges::sort(other_points.begin(), other_points.end(), [&p0](const Point &a, const Point &b) {
    double angle_a = PolarAngle(p0, a);
    double angle_b = PolarAngle(p0, b);
    if (std::abs(angle_a - angle_b) < 1e-10) {
      return DistanceSquared(p0, a) < DistanceSquared(p0, b);
    }
    return angle_a < angle_b;
  });

  return other_points;
}

bool UrinOGrahamPassageSTL::AreAllCollinear(const Point &p0, const std::vector<Point> &points) {
  std::atomic<bool> all_collinear{true};
  unsigned int num_threads = std::thread::hardware_concurrency();
  if (num_threads == 0) {
    num_threads = 2;
  }

  std::vector<std::thread> threads;
  size_t block_size = points.size() / num_threads;

  for (unsigned int thread_idx = 0; thread_idx < num_threads; ++thread_idx) {
    size_t start = (thread_idx == 0) ? 1 : thread_idx * block_size;
    size_t end = (thread_idx == num_threads - 1) ? points.size() : (thread_idx + 1) * block_size;

    threads.emplace_back([&points, &p0, start, end, &all_collinear]() {
      for (size_t i = start; i < end && all_collinear.load(); ++i) {
        if (Orientation(p0, points[0], points[i]) != 0) {
          all_collinear.store(false);
        }
      }
    });
  }

  for (auto &th : threads) {
    th.join();
  }

  return all_collinear.load();
}

std::vector<Point> UrinOGrahamPassageSTL::BuildConvexHull(const Point &p0, const std::vector<Point> &points) {
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

bool UrinOGrahamPassageSTL::RunImpl() {
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

bool UrinOGrahamPassageSTL::PostProcessingImpl() {
  return !GetOutput().empty();
}

}  // namespace urin_o_graham_passage
