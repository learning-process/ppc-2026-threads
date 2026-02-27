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

bool AllPointsSame(const std::vector<Point> &pts) {
  for (size_t i = 1; i < pts.size(); i++) {
    if (pts[i].x != pts[0].x || pts[i].y != pts[0].y) {
      return false;
    }
  }
  return true;
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

  int pivot_idx = 0;
  if (n >= num_threads * 2) {
    std::vector<int> local_pivots(num_threads);
    std::vector<std::thread> threads;
    int chunk = n / num_threads;

    for (int t = 0; t < num_threads; t++) {
      int start = t * chunk;
      int end = (t == num_threads - 1) ? n : (t + 1) * chunk;
      threads.emplace_back([this, &local_pivots, t, start, end]() {
        int local_min = start;
        for (int i = start + 1; i < end; i++) {
          if (points_[i].y < points_[local_min].y ||
              (points_[i].y == points_[local_min].y && points_[i].x < points_[local_min].x)) {
            local_min = i;
          }
        }
        local_pivots[t] = local_min;
      });
    }
    for (auto &th : threads) {
      th.join();
    }

    pivot_idx = local_pivots[0];
    for (int t = 1; t < num_threads; t++) {
      int idx = local_pivots[t];
      if (points_[idx].y < points_[pivot_idx].y ||
          (points_[idx].y == points_[pivot_idx].y && points_[idx].x < points_[pivot_idx].x)) {
        pivot_idx = idx;
      }
    }
  } else {
    for (int i = 1; i < n; i++) {
      if (points_[i].y < points_[pivot_idx].y ||
          (points_[i].y == points_[pivot_idx].y && points_[i].x < points_[pivot_idx].x)) {
        pivot_idx = i;
      }
    }
  }

  std::swap(points_[0], points_[pivot_idx]);
  Point pivot = points_[0];

  auto comparator = [&pivot](const Point &a, const Point &b) {
    double cross = CrossProduct(pivot, a, b);
    if (cross > 0.0) {
      return true;
    }
    if (cross < 0.0) {
      return false;
    }
    return DistSquared(pivot, a) < DistSquared(pivot, b);
  };

  int sort_n = n - 1;

  if (num_threads > 1 && sort_n > num_threads) {
    int chunk = sort_n / num_threads;
    std::vector<std::thread> threads;

    for (int t = 0; t < num_threads; t++) {
      int start = 1 + t * chunk;
      int end = (t == num_threads - 1) ? n : 1 + (t + 1) * chunk;
      threads.emplace_back(
          [this, start, end, &comparator]() { std::sort(points_.begin() + start, points_.begin() + end, comparator); });
    }
    for (auto &th : threads) {
      th.join();
    }

    int merged_end = 1 + chunk;
    for (int t = 1; t < num_threads; t++) {
      int next_end = (t == num_threads - 1) ? n : 1 + (t + 1) * chunk;
      std::inplace_merge(points_.begin() + 1, points_.begin() + merged_end, points_.begin() + next_end, comparator);
      merged_end = next_end;
    }
  } else {
    std::sort(points_.begin() + 1, points_.end(), comparator);
  }

  std::vector<Point> stack;
  stack.reserve(n);
  stack.push_back(points_[0]);
  stack.push_back(points_[1]);

  for (int i = 2; i < n; i++) {
    while (stack.size() > 1 && CrossProduct(stack[stack.size() - 2], stack[stack.size() - 1], points_[i]) <= 0.0) {
      stack.pop_back();
    }
    stack.push_back(points_[i]);
  }

  hull_ = stack;
  return true;
}

bool DergachevAGrahamScanSTL::PostProcessingImpl() {
  GetOutput() = static_cast<int>(hull_.size());
  return true;
}

}  // namespace dergachev_a_graham_scan_stl
