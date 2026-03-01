#include "dergachev_a_graham_scan_stl/stl/include/ops_stl.hpp"

#include <algorithm>
#include <cmath>
#include <thread>
#include <utility>
#include <vector>

#include "dergachev_a_graham_scan_stl/common/include/common.hpp"
#include "util/include/util.hpp"

namespace dergachev_a_graham_scan_stl {

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

bool IsLowerLeft(const Point &a, const Point &b) {
  return a.y < b.y || (a.y == b.y && a.x < b.x);
}

bool AllPointsSame(const std::vector<Point> &pts) {
  int n = static_cast<int>(pts.size());
  for (int i = 1; i < n; i++) {
    if (pts[i].x != pts[0].x || pts[i].y != pts[0].y) {
      return false;
    }
  }
  return true;
}

bool CompareByAngle(const Point &pivot, const Point &a, const Point &b) {
  double cross = CrossProduct(pivot, a, b);
  if (cross > 0.0) {
    return true;
  }
  if (cross < 0.0) {
    return false;
  }
  return DistSquared(pivot, a) < DistSquared(pivot, b);
}

int FindLocalPivot(const std::vector<Point> &pts, int start, int end) {
  int local_min = start;
  for (int i = start + 1; i < end; i++) {
    if (IsLowerLeft(pts[i], pts[local_min])) {
      local_min = i;
    }
  }
  return local_min;
}

int FindPivot(const std::vector<Point> &pts, int num_threads) {
  int n = static_cast<int>(pts.size());
  if (n < num_threads * 2) {
    return FindLocalPivot(pts, 0, n);
  }

  int chunk = n / num_threads;
  std::vector<int> local_pivots(num_threads);
  std::vector<std::thread> threads;

  for (int ti = 0; ti < num_threads; ti++) {
    int start = ti * chunk;
    int end = (ti == num_threads - 1) ? n : ((ti + 1) * chunk);
    threads.emplace_back(
        [&pts, &local_pivots, ti, start, end]() { local_pivots[ti] = FindLocalPivot(pts, start, end); });
  }
  for (auto &th : threads) {
    th.join();
  }

  int pivot_idx = local_pivots[0];
  for (int ti = 1; ti < num_threads; ti++) {
    if (IsLowerLeft(pts[local_pivots[ti]], pts[pivot_idx])) {
      pivot_idx = local_pivots[ti];
    }
  }
  return pivot_idx;
}

void SortByAngle(std::vector<Point> &pts, const Point &pivot, int num_threads) {
  int n = static_cast<int>(pts.size());
  int sort_n = n - 1;

  auto cmp = [&pivot](const Point &a, const Point &b) { return CompareByAngle(pivot, a, b); };

  if (num_threads <= 1 || sort_n <= num_threads) {
    std::sort(pts.begin() + 1, pts.end(), cmp);
    return;
  }

  int chunk = sort_n / num_threads;
  std::vector<std::thread> threads;

  for (int ti = 0; ti < num_threads; ti++) {
    int start = 1 + (ti * chunk);
    int end = (ti == num_threads - 1) ? n : 1 + ((ti + 1) * chunk);
    threads.emplace_back([&pts, start, end, &cmp]() { std::sort(pts.begin() + start, pts.begin() + end, cmp); });
  }
  for (auto &th : threads) {
    th.join();
  }

  int merged_end = 1 + chunk;
  for (int ti = 1; ti < num_threads; ti++) {
    int next_end = (ti == num_threads - 1) ? n : 1 + ((ti + 1) * chunk);
    std::inplace_merge(pts.begin() + 1, pts.begin() + merged_end, pts.begin() + next_end, cmp);
    merged_end = next_end;
  }
}

std::vector<Point> BuildHull(const std::vector<Point> &pts) {
  int n = static_cast<int>(pts.size());
  std::vector<Point> hull;
  hull.reserve(n);
  hull.push_back(pts[0]);
  hull.push_back(pts[1]);

  for (int i = 2; i < n; i++) {
    while (hull.size() > 1 && CrossProduct(hull[hull.size() - 2], hull[hull.size() - 1], pts[i]) <= 0.0) {
      hull.pop_back();
    }
    hull.push_back(pts[i]);
  }
  return hull;
}

}  // namespace

DergachevAGrahamScanSTL::DergachevAGrahamScanSTL(const InType &in) {
  SetTypeOfTask(GetStaticTypeOfTask());
  GetInput() = in;
  GetOutput() = 0;
}

void DergachevAGrahamScanSTL::SetPoints(const std::vector<Point> &pts) {
  points_.assign(pts.begin(), pts.end());
  custom_points_ = true;
}

std::vector<Point> DergachevAGrahamScanSTL::GetHull() const {
  return hull_;
}

bool DergachevAGrahamScanSTL::ValidationImpl() {
  return (GetInput() > 0) && (GetOutput() == 0);
}

bool DergachevAGrahamScanSTL::PreProcessingImpl() {
  int n = GetInput();
  if (!custom_points_) {
    if (n <= 0) {
      points_.clear();
      return true;
    }
    points_.resize(n);
    for (int i = 0; i < n; i++) {
      double angle = 2.0 * kPi * i / n;
      points_[i].x = std::cos(angle);
      points_[i].y = std::sin(angle);
    }
  }
  return true;
}

bool DergachevAGrahamScanSTL::RunImpl() {
  int n = static_cast<int>(points_.size());

  if (n <= 0) {
    hull_.clear();
    return true;
  }
  if (n == 1) {
    hull_ = points_;
    return true;
  }
  if (AllPointsSame(points_)) {
    hull_ = {points_[0]};
    return true;
  }
  if (n == 2) {
    hull_ = points_;
    return true;
  }

  const int num_threads = ppc::util::GetNumThreads();

  int pivot_idx = FindPivot(points_, num_threads);
  std::swap(points_[0], points_[pivot_idx]);

  SortByAngle(points_, points_[0], num_threads);

  hull_ = BuildHull(points_);
  return true;
}

bool DergachevAGrahamScanSTL::PostProcessingImpl() {
  GetOutput() = static_cast<int>(hull_.size());
  return true;
}

}  // namespace dergachev_a_graham_scan_stl
