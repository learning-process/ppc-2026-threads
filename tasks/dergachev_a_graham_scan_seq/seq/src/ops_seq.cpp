#include "dergachev_a_graham_scan_seq/seq/include/ops_seq.hpp"

#include <algorithm>
#include <cmath>
#include <utility>
#include <vector>

#include "dergachev_a_graham_scan_seq/common/include/common.hpp"

namespace dergachev_a_graham_scan_seq {

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

bool AllPointsSame(const std::vector<Point> &pts) {
  for (int i = 1; std::cmp_less(i, pts.size()); i++) {
    if (pts[i].x != pts[0].x || pts[i].y != pts[0].y) {
      return false;
    }
  }
  return true;
}

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

void DergachevAGrahamScanSEQ::SetPoints(const std::vector<Point> &pts) {
  points_.assign(pts.begin(), pts.end());
  custom_points_ = true;
}

std::vector<Point> DergachevAGrahamScanSEQ::GetHull() const {
  return hull_;
}

bool DergachevAGrahamScanSEQ::ValidationImpl() {
  return GetInput() >= 0;
}

bool DergachevAGrahamScanSEQ::PreProcessingImpl() {
  hull_.clear();
  if (!custom_points_) {
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
  }
  return true;
}

bool DergachevAGrahamScanSEQ::RunImpl() {
  hull_.clear();
  std::vector<Point> pts(points_.begin(), points_.end());
  int n = static_cast<int>(pts.size());

  if (n <= 1) {
    hull_ = std::move(pts);
    return true;
  }

  if (AllPointsSame(pts)) {
    hull_.push_back(pts[0]);
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

}  // namespace dergachev_a_graham_scan_seq
