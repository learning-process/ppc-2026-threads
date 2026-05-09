#include "dergachev_a_graham_scan/stl/include/ops_stl.hpp"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <thread>
#include <utility>
#include <vector>

#include "dergachev_a_graham_scan/common/include/common.hpp"
#include "util/include/util.hpp"

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

int FindLocalPivot(const std::vector<Pt> &pts, int start, int end) {
  int best = start;
  for (int i = start + 1; i < end; i++) {
    if (IsLowerLeft(pts[i], pts[best])) {
      best = i;
    }
  }
  return best;
}

int FindPivotParallel(const std::vector<Pt> &pts, int num_threads) {
  int n = static_cast<int>(pts.size());
  if (n < num_threads * 2) {
    return FindLocalPivot(pts, 0, n);
  }

  int chunk = n / num_threads;
  std::vector<int> local_results(num_threads);
  std::vector<std::thread> threads;

  for (int thread_idx = 0; thread_idx < num_threads; thread_idx++) {
    int lo = thread_idx * chunk;
    int hi = (thread_idx == num_threads - 1) ? n : ((thread_idx + 1) * chunk);
    threads.emplace_back(
        [&pts, &local_results, thread_idx, lo, hi]() { local_results[thread_idx] = FindLocalPivot(pts, lo, hi); });
  }
  for (auto &th : threads) {
    th.join();
  }

  int best = local_results[0];
  for (int thread_idx = 1; thread_idx < num_threads; thread_idx++) {
    if (IsLowerLeft(pts[local_results[thread_idx]], pts[best])) {
      best = local_results[thread_idx];
    }
  }
  return best;
}

void GeneratePointsParallel(Pt *data, double step, int n, int num_threads) {
  int chunk = n / num_threads;
  std::vector<std::thread> threads;
  for (int ti = 0; ti < num_threads; ti++) {
    int lo = ti * chunk;
    int hi = (ti == num_threads - 1) ? n : (ti + 1) * chunk;
    threads.emplace_back([data, step, lo, hi]() {
      for (int i = lo; i < hi; i++) {
        data[i] = {std::cos(step * i), std::sin(step * i)};
      }
    });
  }
  for (auto &th : threads) {
    th.join();
  }
}

}  // namespace

DergachevAGrahamScanSTL::DergachevAGrahamScanSTL(const InType &in) {
  SetTypeOfTask(GetStaticTypeOfTask());
  GetInput() = in;
  GetOutput() = 0;
}

bool DergachevAGrahamScanSTL::ValidationImpl() {
  return GetInput() >= 0;
}

bool DergachevAGrahamScanSTL::PreProcessingImpl() {
  hull_.clear();
  int n = GetInput();
  if (n <= 0) {
    points_.clear();
    return true;
  }
  points_.resize(n);
  double step = (2.0 * kPi) / n;

  const int num_threads = ppc::util::GetNumThreads();
  if (num_threads > 1 && n > num_threads * 2) {
    GeneratePointsParallel(points_.data(), step, n, num_threads);
  } else {
    for (int i = 0; i < n; i++) {
      points_[static_cast<std::size_t>(i)] = {std::cos(step * i), std::sin(step * i)};
    }
  }

  if (n > 3) {
    points_.emplace_back(0.0, 0.0);
  }
  return true;
}

bool DergachevAGrahamScanSTL::RunImpl() {
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

  const int num_threads = ppc::util::GetNumThreads();

  int pivot_idx = FindPivotParallel(pts, num_threads);
  std::swap(pts[0], pts[static_cast<std::size_t>(pivot_idx)]);

  Pt pivot = pts[0];
  std::sort(pts.begin() + 1, pts.end(), [&pivot](const Pt &a, const Pt &b) {
    double cross = CrossProduct(pivot, a, b);
    if (cross != 0.0) {
      return cross > 0.0;
    }
    return DistSquared(pivot, a) < DistSquared(pivot, b);
  });

  for (const auto &p : pts) {
    while (hull_.size() > 1 && CrossProduct(hull_[hull_.size() - 2], hull_.back(), p) <= 0.0) {
      hull_.pop_back();
    }
    hull_.push_back(p);
  }

  return true;
}

bool DergachevAGrahamScanSTL::PostProcessingImpl() {
  GetOutput() = static_cast<int>(hull_.size());
  return true;
}

}  // namespace dergachev_a_graham_scan
