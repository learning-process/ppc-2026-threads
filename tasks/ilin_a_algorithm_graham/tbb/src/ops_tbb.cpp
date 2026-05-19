#include "ilin_a_algorithm_graham/tbb/include/ops_tbb.hpp"

#include <tbb/blocked_range.h>
#include <tbb/parallel_reduce.h>
#include <tbb/parallel_sort.h>

#include <cmath>
#include <cstddef>
#include <utility>
#include <vector>

#include "ilin_a_algorithm_graham/common/include/common.hpp"

namespace ilin_a_algorithm_graham {

namespace {
double Orient(const Point &p, const Point &q, const Point &r) {
  return ((q.x - p.x) * (r.y - p.y)) - ((q.y - p.y) * (r.x - p.x));
}

double DistanceSq(const Point &p, const Point &q) {
  double dx = p.x - q.x;
  double dy = p.y - q.y;
  return (dx * dx) + (dy * dy);
}

Point FindLowestLeftmost(const std::vector<Point> &points) {
  return tbb::parallel_reduce(tbb::blocked_range<size_t>(0, points.size()), points[0],
                              [&](const tbb::blocked_range<size_t> &r, Point init) {
    for (size_t i = r.begin(); i < r.end(); ++i) {
      if (points[i].y < init.y || (points[i].y == init.y && points[i].x < init.x)) {
        init = points[i];
      }
    }
    return init;
  }, [](const Point &a, const Point &b) {
    if (a.y < b.y || (a.y == b.y && a.x < b.x)) {
      return a;
    }
    return b;
  });
}

class PointComparator {
 public:
  explicit PointComparator(const Point &p0) : p0_(p0) {}

  bool operator()(const Point &a, const Point &b) const {
    double o = Orient(p0_, a, b);
    if (o != 0.0) {
      return o > 0;
    }
    return DistanceSq(p0_, a) < DistanceSq(p0_, b);
  }

 private:
  Point p0_;
};
}  // namespace

IlinAGrahamTBB::IlinAGrahamTBB(const InType &in) {
  SetTypeOfTask(GetStaticTypeOfTask());
  GetInput() = in;
}

bool IlinAGrahamTBB::ValidationImpl() {
  return !GetInput().points.empty();
}

bool IlinAGrahamTBB::PreProcessingImpl() {
  points_ = GetInput().points;
  hull_.clear();
  return true;
}

bool IlinAGrahamTBB::RunImpl() {
  if (points_.size() < 3) {
    hull_ = points_;
    return true;
  }

  Point p0 = FindLowestLeftmost(points_);

  std::vector<Point> sorted;
  sorted.reserve(points_.size());
  for (const Point &p : points_) {
    if (p.x != p0.x || p.y != p0.y) {
      sorted.push_back(p);
    }
  }

  tbb::parallel_sort(sorted.begin(), sorted.end(), PointComparator(p0));

  std::vector<Point> stack;
  stack.reserve(sorted.size() + 1);
  stack.push_back(p0);
  stack.push_back(sorted[0]);

  for (size_t i = 1; i < sorted.size(); ++i) {
    while (stack.size() >= 2) {
      Point p = stack[stack.size() - 2];
      Point q = stack[stack.size() - 1];
      if (Orient(p, q, sorted[i]) <= 0.0) {
        stack.pop_back();
      } else {
        break;
      }
    }
    stack.push_back(sorted[i]);
  }

  hull_ = std::move(stack);
  return true;
}

bool IlinAGrahamTBB::PostProcessingImpl() {
  GetOutput() = hull_;
  return true;
}

}  // namespace ilin_a_algorithm_graham
