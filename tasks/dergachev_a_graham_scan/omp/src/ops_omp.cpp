#include "dergachev_a_graham_scan/omp/include/ops_omp.hpp"

#include <omp.h>

#include <algorithm>
#include <cmath>
#include <utility>
#include <vector>

#include "dergachev_a_graham_scan/common/include/common.hpp"

namespace dergachev_a_graham_scan {

namespace {

using Pt = std::pair<double, double>;

double CrossProduct(const Pt &o, const Pt &a, const Pt &b) {
  return ((a.first - o.first) * (b.second - o.second)) - ((a.second - o.second) * (b.first - o.first));
}

double DistSquared(const Pt &a, const Pt &b) {
  double dx = a.first - b.first;
  double dy = a.second - b.second;
  return (dx * dx) + (dy * dy);
}

const double kPi = std::acos(-1.0);

bool IsLowerLeft(const Pt &a, const Pt &b) {
  return a.second < b.second || (a.second == b.second && a.first < b.first);
}

int FindPivotIndex(const std::vector<Pt> &pts) {
  int size = static_cast<int>(pts.size());
  int pivot_idx = 0;
#pragma omp parallel default(none) shared(pts, size, pivot_idx)
  {
    int local_idx = 0;
#pragma omp for nowait
    for (int i = 1; i < size; i++) {
      if (IsLowerLeft(pts[i], pts[local_idx])) {
        local_idx = i;
      }
    }
#pragma omp critical
    {
      if (IsLowerLeft(pts[local_idx], pts[pivot_idx])) {
        pivot_idx = local_idx;
      }
    }
  }
  return pivot_idx;
}

void SortByAngle(std::vector<Pt> &pts) {
  Pt pivot = pts[0];
  std::sort(pts.begin() + 1, pts.end(), [&pivot](const Pt &a, const Pt &b) {
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

bool DergachevAGrahamScanOMP::ValidationImpl() {
  return GetInput() >= 0;
}

bool DergachevAGrahamScanOMP::PreProcessingImpl() {
  hull_.clear();
  int n = GetInput();
  if (n <= 0) {
    points_.clear();
    return true;
  }
  points_.resize(n);
  double step = (2.0 * kPi) / n;
  auto *pts_data = points_.data();
#pragma omp parallel for default(none) shared(pts_data, step, n)
  for (int i = 0; i < n; i++) {
    pts_data[i] = {std::cos(step * i), std::sin(step * i)};
  }
  if (n > 3) {
    points_.emplace_back(0.0, 0.0);
  }
  return true;
}

bool DergachevAGrahamScanOMP::RunImpl() {
  hull_.clear();
  std::vector<Pt> pts(points_.begin(), points_.end());
  int n = static_cast<int>(pts.size());

  if (n <= 1 || std::all_of(pts.begin() + 1, pts.end(),
                            [&](const Pt &pt) { return pt.first == pts[0].first && pt.second == pts[0].second; })) {
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

bool DergachevAGrahamScanOMP::PostProcessingImpl() {
  GetOutput() = static_cast<int>(hull_.size());
  return true;
}

}  // namespace dergachev_a_graham_scan
