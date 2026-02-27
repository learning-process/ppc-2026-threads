#include "dergachev_a_graham_scan_tbb/tbb/include/ops_tbb.hpp"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <utility>
#include <vector>

#include "dergachev_a_graham_scan_tbb/common/include/common.hpp"
#include "oneapi/tbb/blocked_range.h"
#include "oneapi/tbb/parallel_for.h"
#include "oneapi/tbb/parallel_reduce.h"
#include "oneapi/tbb/parallel_sort.h"

namespace dergachev_a_graham_scan_tbb {

namespace {

double CrossProduct(const Point &o, const Point &a, const Point &b) {
  return ((a.x - o.x) * (b.y - o.y)) - ((a.y - o.y) * (b.x - o.x));
}

double DistSquared(const Point &a, const Point &b) {
  double dx = a.x - b.x;
  double dy = a.y - b.y;
  return (dx * dx) + (dy * dy);
}

bool IsLower(const Point &a, const Point &b) {
  return a.y < b.y || (a.y == b.y && a.x < b.x);
}

int FindPivot(const std::vector<Point> &pts, int n) {
  return tbb::parallel_reduce(tbb::blocked_range<int>(1, n), 0, [&](const tbb::blocked_range<int> &r, int best) -> int {
    for (int i = r.begin(); i < r.end(); ++i) {
      if (IsLower(pts[static_cast<std::size_t>(i)], pts[static_cast<std::size_t>(best)])) {
        best = i;
      }
    }
    return best;
  }, [&](int a, int b) -> int {
    return IsLower(pts[static_cast<std::size_t>(a)], pts[static_cast<std::size_t>(b)]) ? a : b;
  });
}

bool AllPointsSame(const std::vector<Point> &pts, int n) {
  for (int i = 1; i < n; ++i) {
    if (pts[static_cast<std::size_t>(i)].x != pts[0].x || pts[static_cast<std::size_t>(i)].y != pts[0].y) {
      return false;
    }
  }
  return true;
}

void SortByAngle(std::vector<Point> &pts) {
  Point p0 = pts[0];
  tbb::parallel_sort(pts.begin() + 1, pts.end(), [&p0](const Point &a, const Point &b) {
    double c = CrossProduct(p0, a, b);
    if (c > 0.0) {
      return true;
    }
    if (c < 0.0) {
      return false;
    }
    return DistSquared(p0, a) < DistSquared(p0, b);
  });
}

bool AllCollinear(const std::vector<Point> &pts, int n) {
  for (int i = 2; i < n; ++i) {
    if (CrossProduct(pts[0], pts[1], pts[static_cast<std::size_t>(i)]) != 0.0) {
      return false;
    }
  }
  return true;
}

std::vector<Point> BuildHull(const std::vector<Point> &pts, int n) {
  std::vector<Point> stk;
  stk.reserve(static_cast<std::size_t>(n));
  stk.push_back(pts[0]);
  stk.push_back(pts[1]);

  for (int i = 2; i < n; ++i) {
    while (stk.size() >= 2 &&
           CrossProduct(stk[stk.size() - 2], stk[stk.size() - 1], pts[static_cast<std::size_t>(i)]) <= 0.0) {
      stk.pop_back();
    }
    stk.push_back(pts[static_cast<std::size_t>(i)]);
  }
  return stk;
}

}  // namespace

DergachevAGrahamScanTBB::DergachevAGrahamScanTBB(const InType &in) {
  SetTypeOfTask(GetStaticTypeOfTask());
  GetInput() = in;
  GetOutput() = 0;
}

void DergachevAGrahamScanTBB::SetPoints(const std::vector<Point> &pts) {
  points_ = pts;
  custom_points_ = true;
}

std::vector<Point> DergachevAGrahamScanTBB::GetHull() const {
  return hull_;
}

bool DergachevAGrahamScanTBB::ValidationImpl() {
  return GetInput() >= 0;
}

bool DergachevAGrahamScanTBB::PreProcessingImpl() {
  hull_.clear();
  GetOutput() = 0;

  if (!custom_points_) {
    int n = GetInput();
    if (n <= 0) {
      points_.clear();
      return true;
    }
    points_.resize(static_cast<std::size_t>(n));
    const double pi2 = 2.0 * std::acos(-1.0);
    auto dn = static_cast<double>(n);
    tbb::parallel_for(0, n, [&](int i) {
      double angle = (pi2 * static_cast<double>(i)) / dn;
      points_[static_cast<std::size_t>(i)] = {.x = std::cos(angle), .y = std::sin(angle)};
    });
  }

  return true;
}

bool DergachevAGrahamScanTBB::RunImpl() {
  int n = static_cast<int>(points_.size());

  if (n <= 0) {
    hull_.clear();
    GetOutput() = 0;
    return true;
  }

  if (n == 1) {
    hull_ = points_;
    GetOutput() = 1;
    return true;
  }

  int pivot = FindPivot(points_, n);
  if (pivot != 0) {
    std::swap(points_[0], points_[static_cast<std::size_t>(pivot)]);
  }

  if (AllPointsSame(points_, n)) {
    hull_ = {points_[0]};
    GetOutput() = 1;
    return true;
  }

  SortByAngle(points_);

  if (AllCollinear(points_, n)) {
    hull_ = {points_[0], points_[static_cast<std::size_t>(n - 1)]};
    GetOutput() = 2;
    return true;
  }

  hull_ = BuildHull(points_, n);
  GetOutput() = static_cast<OutType>(hull_.size());
  return true;
}

bool DergachevAGrahamScanTBB::PostProcessingImpl() {
  return std::cmp_equal(GetOutput(), hull_.size());
}

}  // namespace dergachev_a_graham_scan_tbb
