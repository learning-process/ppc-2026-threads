#include "urin_o_graham_passage/all/include/ops_all.hpp"

#include <algorithm>
#include <atomic>
#include <cmath>
#include <cstddef>
// #include <thread>
#include <vector>

// Подключение TBB
#include <tbb/blocked_range.h>
#include <tbb/concurrent_vector.h>
#include <tbb/parallel_for.h>
#include <tbb/parallel_reduce.h>

// OpenMP (используется для демонстрации в одном из этапов)
#include <omp.h>

#include "urin_o_graham_passage/common/include/common.hpp"
// #include "util/include/util.hpp"

namespace urin_o_graham_passage {

UrinOGrahamPassageALL::UrinOGrahamPassageALL(const InType &in) {
  SetTypeOfTask(GetStaticTypeOfTask());
  GetInput() = in;
  GetOutput() = OutType();
}

bool UrinOGrahamPassageALL::ValidationImpl() {
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

bool UrinOGrahamPassageALL::PreProcessingImpl() {
  GetOutput().clear();
  return true;
}

Point UrinOGrahamPassageALL::FindLowestPoint(const InType &points) {
  Point lowest = points[0];
  for (size_t i = 1; i < points.size(); ++i) {
    if (points[i].y < lowest.y - 1e-10 ||
        (std::abs(points[i].y - lowest.y) < 1e-10 && points[i].x < lowest.x - 1e-10)) {
      lowest = points[i];
    }
  }
  return lowest;
}

Point UrinOGrahamPassageALL::FindLowestPointParallel(const InType &points) {
  return tbb::parallel_reduce(tbb::blocked_range<size_t>(0, points.size()), points[0],
                              [&points](const tbb::blocked_range<size_t> &range, Point current_min) {
    for (size_t i = range.begin(); i < range.end(); ++i) {
      if (points[i].y < current_min.y - 1e-10 ||
          (std::abs(points[i].y - current_min.y) < 1e-10 && points[i].x < current_min.x - 1e-10)) {
        current_min = points[i];
      }
    }
    return current_min;
  }, [](const Point &a, const Point &b) {
    if (a.y < b.y - 1e-10 || (std::abs(a.y - b.y) < 1e-10 && a.x < b.x - 1e-10)) {
      return a;
    }
    return b;
  });
}

double UrinOGrahamPassageALL::PolarAngle(const Point &base, const Point &p) {
  double dx = p.x - base.x;
  double dy = p.y - base.y;
  if (std::abs(dx) < 1e-10 && std::abs(dy) < 1e-10) {
    return -1e10;
  }
  return std::atan2(dy, dx);
}

int UrinOGrahamPassageALL::Orientation(const Point &p, const Point &q, const Point &r) {
  double val = ((q.x - p.x) * (r.y - p.y)) - ((q.y - p.y) * (r.x - p.x));
  if (std::abs(val) < 1e-10) {
    return 0;
  }
  return (val > 0) ? 1 : -1;
}

double UrinOGrahamPassageALL::DistanceSquared(const Point &p1, const Point &p2) {
  double dx = p2.x - p1.x;
  double dy = p2.y - p1.y;
  return (dx * dx) + (dy * dy);
}

std::vector<Point> UrinOGrahamPassageALL::PrepareOtherPointsParallel(const InType &points, const Point &p0) {
  tbb::concurrent_vector<Point> other_points_concurrent;
  tbb::parallel_for(tbb::blocked_range<size_t>(0, points.size()),
                    [&points, &p0, &other_points_concurrent](const tbb::blocked_range<size_t> &range) {
    for (size_t i = range.begin(); i < range.end(); ++i) {
      if (points[i] != p0) {
        other_points_concurrent.push_back(points[i]);
      }
    }
  });
  std::vector<Point> other_points(other_points_concurrent.begin(), other_points_concurrent.end());
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

bool UrinOGrahamPassageALL::AreAllCollinear(const Point &p0, const std::vector<Point> &points) {
  std::atomic<bool> all_collinear{true};
  tbb::parallel_for(tbb::blocked_range<size_t>(1, points.size()),
                    [&points, &p0, &all_collinear](const tbb::blocked_range<size_t> &range) {
    for (size_t i = range.begin(); i < range.end() && all_collinear.load(); ++i) {
      if (Orientation(p0, points[0], points[i]) != 0) {
        all_collinear.store(false);
      }
    }
  });
  return all_collinear.load();
}

std::vector<Point> UrinOGrahamPassageALL::BuildConvexHull(const Point &p0, const std::vector<Point> &points) {
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

bool UrinOGrahamPassageALL::RunImpl() {
  const InType &points = GetInput();
  if (points.size() < 3) {
    return false;
  }

  // Поиск самой нижней точки – через TBB parallel_reduce
  Point p0 = FindLowestPointParallel(points);

  // Сбор точек, отличных от p0 – через TBB concurrent_vector
  std::vector<Point> other_points = PrepareOtherPointsParallel(points, p0);
  if (other_points.empty()) {
    return false;
  }

  // Проверка на коллинеарность – через TBB parallel_for с атомарным флагом
  if (AreAllCollinear(p0, other_points)) {
    GetOutput() = {p0, other_points.back()};
    return true;
  }

  // Построение выпуклой оболочки – последовательно
  GetOutput() = BuildConvexHull(p0, other_points);
  return true;
}

bool UrinOGrahamPassageALL::PostProcessingImpl() {
  return !GetOutput().empty();
}

}  // namespace urin_o_graham_passage
