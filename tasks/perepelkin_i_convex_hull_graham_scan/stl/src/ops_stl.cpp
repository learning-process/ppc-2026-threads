#include "perepelkin_i_convex_hull_graham_scan/stl/include/ops_stl.hpp"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <thread>
#include <utility>
#include <vector>

#include "perepelkin_i_convex_hull_graham_scan/common/include/common.hpp"
#include "util/include/util.hpp"

namespace perepelkin_i_convex_hull_graham_scan {

PerepelkinIConvexHullGrahamScanSTL::PerepelkinIConvexHullGrahamScanSTL(const InType &in) {
  SetTypeOfTask(GetStaticTypeOfTask());
  GetInput() = in;
  GetOutput() = std::vector<std::pair<double, double>>();
}

bool PerepelkinIConvexHullGrahamScanSTL::ValidationImpl() {
  return GetOutput().empty();
}

bool PerepelkinIConvexHullGrahamScanSTL::PreProcessingImpl() {
  return true;
}

bool PerepelkinIConvexHullGrahamScanSTL::RunImpl() {
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
  HullConstruction(hull, pts, pivot);

  GetOutput() = std::move(hull);
  return true;
}

size_t PerepelkinIConvexHullGrahamScanSTL::FindPivotParallel(const std::vector<std::pair<double, double>> &pts) {
  const int threads = std::min(ppc::util::GetNumThreads(), static_cast<int>(pts.size()));

  // Partitioning
  std::vector<int> start(threads + 1);
  DataPartitioning(pts.size(), threads, start);

  // Parallel search
  std::vector<size_t> local_idx(threads, 0);
  {
    std::vector<std::jthread> workers(threads);

    for (int tid = 0; tid < threads; tid++) {
      workers.emplace_back([&, tid]() {
        size_t begin = start[tid];
        size_t end = start[tid + 1];

        size_t local = begin;

        for (size_t i = begin + 1; i < end; i++) {
          if (pts[i].second < pts[local].second ||
              (pts[i].second == pts[local].second && pts[i].first < pts[local].first)) {
            local = i;
          }
        }

        local_idx[tid] = local;
      });
    }
  }

  // Sequential reduction
  size_t pivot_idx = 0;

  for (int tid = 0; tid < threads; tid++) {
    size_t i = local_idx[tid];

    if (pts[i].second < pts[pivot_idx].second ||
        (pts[i].second == pts[pivot_idx].second && pts[i].first < pts[pivot_idx].first)) {
      pivot_idx = i;
    }
  }

  return pivot_idx;
}

void PerepelkinIConvexHullGrahamScanSTL::ParallelSort(std::vector<std::pair<double, double>> &data,
                                                      const std::pair<double, double> &pivot) {
  const int threads = std::min(ppc::util::GetNumThreads(), static_cast<int>(data.size()));

  // Partitioning
  std::vector<int> start(threads + 1);
  DataPartitioning(data.size(), threads, start);

  // Parallel local sorting
  {
    std::vector<std::jthread> workers(threads);

    for (int tid = 0; tid < threads; tid++) {
      workers.emplace_back([&, tid]() {
        std::sort(data.begin() + start[tid], data.begin() + start[tid + 1],
                  [&](const auto &a, const auto &b) { return AngleCmp(a, b, pivot); });
      });
    }
  }

  // Merge sorted segments
  for (int size = 1; size < threads; size *= 2) {
    std::vector<std::jthread> workers(threads);

    for (int i = 0; i < threads; i += 2 * size) {
      if (i + size >= threads) {
        continue;
      }

      workers.emplace_back([&, i, size]() {
        int left = start[i];
        int mid = start[i + size];
        int right = start[std::min(i + (2 * size), threads)];

        std::inplace_merge(data.begin() + left, data.begin() + mid, data.begin() + right,
                           [&](const auto &a, const auto &b) { return AngleCmp(a, b, pivot); });
      });
    }
  }
}

void PerepelkinIConvexHullGrahamScanSTL::DataPartitioning(size_t total_size, const int &threads,
                                                          std::vector<int> &start) {
  size_t base = total_size / threads;
  size_t rem = total_size % threads;

  size_t offset = 0;

  for (int i = 0; i < threads; i++) {
    start[i] = static_cast<int>(offset);

    size_t extra = std::cmp_less(i, rem) ? 1 : 0;
    offset += base + extra;
  }

  start[threads] = static_cast<int>(total_size);
}

void PerepelkinIConvexHullGrahamScanSTL::HullConstruction(std::vector<std::pair<double, double>> &hull,
                                                          const std::vector<std::pair<double, double>> &pts,
                                                          const std::pair<double, double> &pivot) {
  hull.reserve(pts.size() + 1);

  hull.push_back(pivot);
  hull.push_back(pts[0]);

  for (size_t i = 1; i < pts.size(); i++) {
    while (hull.size() >= 2 && Orientation(hull[hull.size() - 2], hull[hull.size() - 1], pts[i]) <= 0) {
      hull.pop_back();
    }

    hull.push_back(pts[i]);
  }
}

double PerepelkinIConvexHullGrahamScanSTL::Orientation(const std::pair<double, double> &p,
                                                       const std::pair<double, double> &q,
                                                       const std::pair<double, double> &r) {
  double val = ((q.first - p.first) * (r.second - p.second)) - ((q.second - p.second) * (r.first - p.first));

  if (std::abs(val) < 1e-9) {
    return 0.0;
  }

  return val;
}

bool PerepelkinIConvexHullGrahamScanSTL::AngleCmp(const std::pair<double, double> &a,
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

bool PerepelkinIConvexHullGrahamScanSTL::PostProcessingImpl() {
  return true;
}

}  // namespace perepelkin_i_convex_hull_graham_scan
