#include "peterson_r_graham_scan_omp/omp/include/ops_omp.hpp"

#include <omp.h>

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <numbers>
#include <vector>

#include "peterson_r_graham_scan_omp/common/include/common.hpp"

namespace peterson_r_graham_scan_omp {

namespace {

double CalculateDeterminant(const Point &origin, const Point &a, const Point &b) {
  return ((a.x - origin.x) * (b.y - origin.y)) - ((a.y - origin.y) * (b.x - origin.x));
}

double GetSquaredDist(const Point &a, const Point &b) {
  const double dx = a.x - b.x;
  const double dy = a.y - b.y;
  return (dx * dx) + (dy * dy);
}

bool AreAllPointsIdentical(const std::vector<Point> &pts) {
  if (pts.empty()) {
    return true;
  }
  bool identical = true;
  const int n = static_cast<int>(pts.size());
#pragma omp parallel for default(none) shared(pts, n) reduction(&& : identical)
  for (int i = 1; i < n; ++i) {
    if (pts[i].x != pts[0].x || pts[i].y != pts[0].y) {
      identical = false;
    }
  }
  return identical;
}

std::size_t LocatePivotIndex(const std::vector<Point> &pts) {
  const int n = static_cast<int>(pts.size());
  int best_global_idx = 0;

#pragma omp parallel default(none) shared(pts, n, best_global_idx)
  {
    int thread_best_idx = 0;
#pragma omp for nowait
    for (int i = 1; i < n; ++i) {
      if (pts[i].y < pts[thread_best_idx].y ||
          (pts[i].y == pts[thread_best_idx].y && pts[i].x < pts[thread_best_idx].x)) {
        thread_best_idx = i;
      }
    }
#pragma omp critical
    {
      if (pts[thread_best_idx].y < pts[best_global_idx].y ||
          (pts[thread_best_idx].y == pts[best_global_idx].y && pts[thread_best_idx].x < pts[best_global_idx].x)) {
        best_global_idx = thread_best_idx;
      }
    }
  }
  return static_cast<std::size_t>(best_global_idx);
}

void SortPointsByPolarAngle(std::vector<Point> &pts) {
  Point pivot = pts[0];
  std::sort(pts.begin() + 1, pts.end(), [pivot](const Point &lhs, const Point &rhs) {
    const double det = CalculateDeterminant(pivot, lhs, rhs);
    if (det > 0.0) {
      return true;
    }
    if (det < 0.0) {
      return false;
    }
    return GetSquaredDist(pivot, lhs) < GetSquaredDist(pivot, rhs);
  });
}

}  // namespace

PetersonRGrahamScanOMP::PetersonRGrahamScanOMP(const InType &in) {
  SetTypeOfTask(GetStaticTypeOfTask());
  GetInput() = in;
  GetOutput() = 0;
}

void PetersonRGrahamScanOMP::SetPoints(const std::vector<Point> &pts) {
  input_points_ = pts;
  is_custom_input_ = true;
}

std::vector<Point> PetersonRGrahamScanOMP::GetHull() const {
  return convex_hull_;
}

bool PetersonRGrahamScanOMP::ValidationImpl() {
  return GetInput() >= 0;
}

bool PetersonRGrahamScanOMP::PreProcessingImpl() {
  convex_hull_.clear();

  if (!is_custom_input_) {
    input_points_.clear();

    if (GetInput() <= 0) {
      return true;
    }

    const int point_count = GetInput();
    input_points_.reserve(point_count);
    for (int i = 0; i < point_count; ++i) {
      const double angle = (2.0 * std::numbers::pi * static_cast<double>(i)) / static_cast<double>(point_count);
      input_points_.emplace_back(Point{std::cos(angle), std::sin(angle)});
    }
  } else {
    if (input_points_.size() != static_cast<std::size_t>(GetInput())) {
      GetInput() = static_cast<InType>(input_points_.size());
    }
  }

  return true;
}

bool PetersonRGrahamScanOMP::RunImpl() {
  convex_hull_.clear();

  const int n = static_cast<int>(input_points_.size());
  if (n == 0) {
    return true;
  }

  if (AreAllPointsIdentical(input_points_)) {
    convex_hull_.push_back(input_points_[0]);
    return true;
  }

  if (n < 3) {
    convex_hull_ = input_points_;
    return true;
  }

  const std::size_t pivot_idx = LocatePivotIndex(input_points_);
  std::iter_swap(input_points_.begin(), input_points_.begin() + pivot_idx);
  SortPointsByPolarAngle(input_points_);

  convex_hull_.push_back(input_points_[0]);
  convex_hull_.push_back(input_points_[1]);

  for (int i = 2; i < n; ++i) {
    while (convex_hull_.size() >= 2) {
      const Point &p1 = convex_hull_[convex_hull_.size() - 2];
      const Point &p2 = convex_hull_.back();
      if (CalculateDeterminant(p1, p2, input_points_[i]) <= 0.0) {
        convex_hull_.pop_back();
      } else {
        break;
      }
    }
    convex_hull_.push_back(input_points_[i]);
  }

  return true;
}

bool PetersonRGrahamScanOMP::PostProcessingImpl() {
  GetOutput() = static_cast<OutType>(convex_hull_.size());
  return true;
}

}  // namespace peterson_r_graham_scan_omp
