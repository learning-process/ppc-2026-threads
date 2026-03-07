#include "dergachev_a_graham_scan/seq/include/ops_seq.hpp"

#include <algorithm>
#include <cmath>
#include <utility>
#include <vector>

#include "dergachev_a_graham_scan/common/include/common.hpp"

namespace dergachev_a_graham_scan {

namespace {

double CrossProduct(const Point &o, const Point &a, const Point &b) {
  return ((a.x - o.x) * (b.y - o.y)) - ((a.y - o.y) * (b.x - o.x));
}

double DistSquared(const Point &a, const Point &b) {
  double dx = a.x - b.x;
  double dy = a.y - b.y;
  return (dx * dx) + (dy * dy);
}

const double kPi = std::acos(-1.0);

int FindPivotIndex(const std::vector<Point> &pts) {
  int pivot_idx = 0;
  for (int i = 1; std::cmp_less(i, pts.size()); i++) {
    if (pts[i].y < pts[pivot_idx].y || (pts[i].y == pts[pivot_idx].y && pts[i].x < pts[pivot_idx].x)) {
      pivot_idx = i;
    }
  }
  return pivot_idx;
}

void SortByAngle(std::vector<Point> &pts) {
  Point pivot = pts[0];
  std::sort(pts.begin() + 1, pts.end(), [&pivot](const Point &a, const Point &b) {
    double cross = CrossProduct(pivot, a, b);
    if (cross > 0.0) {
      return true;
    }
    if (cross < 0.0) {
      return false;
    }
    return DistSquared(pivot, a) < DistSquared(pivot, b);
  });
}

}  // namespace

DergachevAGrahamScanSEQ::DergachevAGrahamScanSEQ(const InType &in) {
  SetTypeOfTask(GetStaticTypeOfTask());
  GetInput() = in;
  GetOutput() = 0;
}

bool DergachevAGrahamScanSEQ::ValidationImpl() {
  return GetInput() >= 0;
}

bool DergachevAGrahamScanSEQ::PreProcessingImpl() {
  hull_.clear();
  int n = GetInput();
  if (n <= 0) {
    points_.clear();
    return true;
  }
  points_.resize(n);
  double step = (2.0 * kPi) / n;
  for (int i = 0; i < n; i++) {
    points_[i] = {.x = std::cos(step * i), .y = std::sin(step * i)};
  }
  if (n > 3) {
    points_.push_back({.x = 0.0, .y = 0.0});
  }
  return true;
}

bool DergachevAGrahamScanSEQ::RunImpl() {
  hull_.clear();
  std::vector<Point> pts(points_.begin(), points_.end());
  int n = static_cast<int>(pts.size());

  if (n <= 1 ||
      std::all_of(pts.begin() + 1, pts.end(), [&](const Point &pt) { return pt.x == pts[0].x && pt.y == pts[0].y; })) {
    if (!pts.empty()) {
      hull_.push_back(pts[0]);
    }
    return true;
  }

  int pivot_idx = FindPivotIndex(pts);
  std::swap(pts[0], pts[pivot_idx]);
  SortByAngle(pts);

  for (const auto &p : pts) {
    while (hull_.size() > 1 && CrossProduct(hull_[hull_.size() - 2], hull_.back(), p) <= 0.0) {
      hull_.pop_back();
    }
    hull_.push_back(p);
  }

  return true;
}

bool DergachevAGrahamScanSEQ::PostProcessingImpl() {
  GetOutput() = static_cast<int>(hull_.size());
  return true;
}

}  // namespace dergachev_a_graham_scan
