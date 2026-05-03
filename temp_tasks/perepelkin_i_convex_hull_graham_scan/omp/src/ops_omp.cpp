#include "perepelkin_i_convex_hull_graham_scan/omp/include/ops_omp.hpp"

#include <omp.h>

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <utility>
#include <vector>

#include "perepelkin_i_convex_hull_graham_scan/common/include/common.hpp"
#include "util/include/util.hpp"

namespace perepelkin_i_convex_hull_graham_scan {

PerepelkinIConvexHullGrahamScanOMP::PerepelkinIConvexHullGrahamScanOMP(const InType &in) {
  SetTypeOfTask(GetStaticTypeOfTask());
  GetInput() = in;
  GetOutput() = std::vector<std::pair<double, double>>();
}

bool PerepelkinIConvexHullGrahamScanOMP::ValidationImpl() {
  return GetOutput().empty();
}

bool PerepelkinIConvexHullGrahamScanOMP::PreProcessingImpl() {
  return true;
}
bool PerepelkinIConvexHullGrahamScanOMP::RunImpl() {
  const auto &data = GetInput();

  if (data.size() < 2) {
    GetOutput() = data;
    return true;
  }

  std::vector<std::pair<double, double>> pts = data;

  // Find pivot
  size_t pivot_idx = FindPivotParallel(pts);

  std::pair<double, double> pivot = pts[pivot_idx];
  pts.erase(pts.begin() + static_cast<std::ptrdiff_t>(pivot_idx));

  // Parallel sorting
  ParallelSort(pts, pivot);

  // Sequential hull construction
  std::vector<std::pair<double, double>> hull;
  hull.reserve(pts.size() + 1);

  hull.push_back(pivot);
  hull.push_back(pts[0]);

  for (size_t i = 1; i < pts.size(); i++) {
    while (hull.size() >= 2 && Orientation(hull[hull.size() - 2], hull[hull.size() - 1], pts[i]) <= 0) {
      hull.pop_back();
    }

    hull.push_back(pts[i]);
  }

  GetOutput() = std::move(hull);
  return true;
}

size_t PerepelkinIConvexHullGrahamScanOMP::FindPivotParallel(const std::vector<std::pair<double, double>> &pts) {
  size_t pivot_idx = 0;

#pragma omp parallel default(none) shared(pts, pivot_idx) num_threads(ppc::util::GetNumThreads())
  {
    size_t local_idx = pivot_idx;

#pragma omp for nowait
    for (size_t i = 1; i < pts.size(); i++) {
      if (pts[i].second < pts[local_idx].second ||
          (pts[i].second == pts[local_idx].second && pts[i].first < pts[local_idx].first)) {
        local_idx = i;
      }
    }

#pragma omp critical
    {
      if (pts[local_idx].second < pts[pivot_idx].second ||
          (pts[local_idx].second == pts[pivot_idx].second && pts[local_idx].first < pts[pivot_idx].first)) {
        pivot_idx = local_idx;
      }
    }
  }

  return pivot_idx;
}

void PerepelkinIConvexHullGrahamScanOMP::ParallelSort(std::vector<std::pair<double, double>> &data,
                                                      const std::pair<double, double> &pivot) {
  size_t n = data.size();
  int threads = ppc::util::GetNumThreads();

  if (n < 10000) {
    std::ranges::sort(data, [&](const auto &a, const auto &b) { return AngleCmp(a, b, pivot); });
    return;
  }

  std::vector<int> start(threads + 1);
  for (int i = 0; i <= threads; i++) {
    start[i] = static_cast<int>(i * n / threads);
  }

#pragma omp parallel default(none) shared(data, start, pivot) num_threads(ppc::util::GetNumThreads())
  {
    int tid = omp_get_thread_num();
    std::sort(data.begin() + start[tid], data.begin() + start[tid + 1],
              [&](const auto &a, const auto &b) { return AngleCmp(a, b, pivot); });
  }

  // Merge sorted segments
  for (int size = 1; size < threads; size *= 2) {
#pragma omp parallel for default(none) shared(data, start, threads, size, pivot) num_threads(ppc::util::GetNumThreads())
    for (int i = 0; i < threads; i += 2 * size) {
      if (i + size >= threads) {
        continue;
      }

      int left = start[i];
      int mid = start[i + size];
      int right = start[std::min(i + (2 * size), threads)];

      std::inplace_merge(data.begin() + left, data.begin() + mid, data.begin() + right,
                         [&](const auto &a, const auto &b) { return AngleCmp(a, b, pivot); });
    }
  }
}

double PerepelkinIConvexHullGrahamScanOMP::Orientation(const std::pair<double, double> &p,
                                                       const std::pair<double, double> &q,
                                                       const std::pair<double, double> &r) {
  double val = ((q.first - p.first) * (r.second - p.second)) - ((q.second - p.second) * (r.first - p.first));

  if (std::abs(val) < 1e-9) {
    return 0.0;
  }

  return val;
}

bool PerepelkinIConvexHullGrahamScanOMP::AngleCmp(const std::pair<double, double> &a,
                                                  const std::pair<double, double> &b,
                                                  const std::pair<double, double> &pivot) {
  double dx1 = a.first - pivot.first;
  double dy1 = a.second - pivot.second;
  double dx2 = b.first - pivot.first;
  double dy2 = b.second - pivot.second;

  double cross = (dx1 * dy2) - (dy1 * dx2);

  if (std::abs(cross) < 1e-9) {
    double dist1 = (dx1 * dx1) + (dy1 * dy1);
    double dist2 = (dx2 * dx2) + (dy2 * dy2);
    return dist1 < dist2;
  }

  return cross > 0;
}

bool PerepelkinIConvexHullGrahamScanOMP::PostProcessingImpl() {
  return true;
}

}  // namespace perepelkin_i_convex_hull_graham_scan
