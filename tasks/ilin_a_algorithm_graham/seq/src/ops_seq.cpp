#include "ilin_a_algorithm_graham/seq/include/ops_seq.hpp"

#include <algorithm>
#include <cmath>
#include <vector>

namespace ilin_a_algorithm_graham {

namespace {
double Orient(const Point& p, const Point& q, const Point& r) {
  return (q.x - p.x) * (r.y - p.y) - (q.y - p.y) * (r.x - p.x);
}

double DistanceSq(const Point& p, const Point& q) {
  double dx = p.x - q.x;
  double dy = p.y - q.y;
  return dx * dx + dy * dy;
}

Point FindLowestLeftmost(const std::vector<Point>& points) {
  Point p0 = points[0];
  for (size_t i = 1; i < points.size(); ++i) {
    if (points[i].y < p0.y || (points[i].y == p0.y && points[i].x < p0.x)) {
      p0 = points[i];
    }
  }
  return p0;
}

class PointComparator {
 public:
  explicit PointComparator(const Point& p0) : p0_(p0) {}
  
  bool operator()(const Point& a, const Point& b) const {
    double o = Orient(p0_, a, b);
    if (o != 0.0) return o > 0;
    return DistanceSq(p0_, a) < DistanceSq(p0_, b);
  }

 private:
  Point p0_;
};
}  // namespace

IlinAGrahamSEQ::IlinAGrahamSEQ(const InType &in) {
  SetTypeOfTask(GetStaticTypeOfTask());
  GetInput() = in;
}

bool IlinAGrahamSEQ::ValidationImpl() {
  return !GetInput().points.empty();
}

bool IlinAGrahamSEQ::PreProcessingImpl() {
  points_ = GetInput().points;
  hull_.clear();
  return true;
}

bool IlinAGrahamSEQ::RunImpl() {
  if (points_.size() < 3) {
    hull_ = points_;
    return true;
  }

  Point p0 = FindLowestLeftmost(points_);
  
  std::vector<Point> sorted;
  for (const Point& p : points_) {
    if (p.x != p0.x || p.y != p0.y) {
      sorted.push_back(p);
    }
  }
  
  std::sort(sorted.begin(), sorted.end(), PointComparator(p0));

  std::vector<Point> stack;
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

  hull_ = stack;
  return true;
}

bool IlinAGrahamSEQ::PostProcessingImpl() {
  GetOutput() = hull_;
  return true;
}

}  // namespace ilin_a_algorithm_graham