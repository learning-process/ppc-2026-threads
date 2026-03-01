#include "dergachev_a_graham_scan_omp/omp/include/ops_omp.hpp"

#include <omp.h>

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <vector>

#include "dergachev_a_graham_scan_omp/common/include/common.hpp"

namespace dergachev_a_graham_scan_omp {

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
  if (pts.empty()) {
    return true;
  }
  bool all_same = true;
  int size = static_cast<int>(pts.size());
#pragma omp parallel for default(none) shared(pts, size) reduction(&& : all_same)
  for (int i = 1; i < size; i++) {
    if (pts[i].x != pts[0].x || pts[i].y != pts[0].y) {
      all_same = false;
    }
  }
  return all_same;
}

std::size_t FindPivotIndex(const std::vector<Point> &pts) {
  int size = static_cast<int>(pts.size());
  int pivot_idx = 0;
#pragma omp parallel default(none) shared(pts, size, pivot_idx)
  {
    int local_idx = 0;
#pragma omp for nowait
    for (int i = 1; i < size; i++) {
      if (pts[i].y < pts[local_idx].y || (pts[i].y == pts[local_idx].y && pts[i].x < pts[local_idx].x)) {
        local_idx = i;
      }
    }
#pragma omp critical
    {
      if (pts[local_idx].y < pts[pivot_idx].y ||
          (pts[local_idx].y == pts[pivot_idx].y && pts[local_idx].x < pts[pivot_idx].x)) {
        pivot_idx = local_idx;
      }
    }
  }
  return static_cast<std::size_t>(pivot_idx);
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

DergachevAGrahamScanOMP::DergachevAGrahamScanOMP(const InType &in) {
  SetTypeOfTask(GetStaticTypeOfTask());
  GetInput() = in;
  GetOutput() = 0;
}

void DergachevAGrahamScanOMP::SetPoints(const std::vector<Point> &pts) {
  points_ = pts;
  custom_points_ = true;
}

std::vector<Point> DergachevAGrahamScanOMP::GetHull() const {
  return hull_;
}

bool DergachevAGrahamScanOMP::ValidationImpl() {
  return GetInput() >= 0;
}

bool DergachevAGrahamScanOMP::PreProcessingImpl() {
  hull_.clear();

  if (!custom_points_) {
    points_.clear();

    if (GetInput() <= 0) {
      return true;
    }

    const int n = GetInput();
    points_.reserve(n);
    for (int i = 0; i < n; ++i) {
      double angle = (2.0 * kPi * static_cast<double>(i)) / static_cast<double>(n);
      Point p{.x = std::cos(angle), .y = std::sin(angle)};
      points_.push_back(p);
    }
  } else {
    if (points_.size() != static_cast<std::size_t>(GetInput())) {
      GetInput() = static_cast<InType>(points_.size());
    }
  }

  return true;
}

bool DergachevAGrahamScanOMP::RunImpl() {
  hull_.clear();

  const int n = static_cast<int>(points_.size());
  if (n == 0) {
    return true;
  }

  if (AllPointsSame(points_)) {
    hull_.push_back(points_[0]);
    return true;
  }

  if (n == 1) {
    hull_.push_back(points_[0]);
    return true;
  }

  if (n == 2) {
    hull_ = points_;
    return true;
  }

  std::size_t pivot_idx = FindPivotIndex(points_);
  std::swap(points_[0], points_[pivot_idx]);
  SortByAngle(points_);

  hull_.push_back(points_[0]);
  hull_.push_back(points_[1]);

  for (int i = 2; i < n; ++i) {
    while (hull_.size() >= 2) {
      const Point &p1 = hull_[hull_.size() - 2];
      const Point &p2 = hull_[hull_.size() - 1];
      if (CrossProduct(p1, p2, points_[i]) <= 0.0) {
        hull_.pop_back();
      } else {
        break;
      }
    }
    hull_.push_back(points_[i]);
  }

  return true;
}

bool DergachevAGrahamScanOMP::PostProcessingImpl() {
  GetOutput() = static_cast<OutType>(hull_.size());
  return true;
}

}  // namespace dergachev_a_graham_scan_omp
