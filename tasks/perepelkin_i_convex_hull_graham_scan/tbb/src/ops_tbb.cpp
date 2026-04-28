#include "perepelkin_i_convex_hull_graham_scan/tbb/include/ops_tbb.hpp"

#include <tbb/blocked_range.h>
#include <tbb/global_control.h>
#include <tbb/parallel_reduce.h>
#include <tbb/parallel_sort.h>

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <utility>
#include <vector>

#include "perepelkin_i_convex_hull_graham_scan/common/include/common.hpp"
#include "util/include/util.hpp"

namespace perepelkin_i_convex_hull_graham_scan {

PerepelkinIConvexHullGrahamScanTBB::PerepelkinIConvexHullGrahamScanTBB(const InType &in) {
  SetTypeOfTask(GetStaticTypeOfTask());
  GetInput() = in;
  GetOutput() = std::vector<std::pair<double, double>>();
}

bool PerepelkinIConvexHullGrahamScanTBB::ValidationImpl() {
  return GetOutput().empty();
}

bool PerepelkinIConvexHullGrahamScanTBB::PreProcessingImpl() {
  return true;
}

bool PerepelkinIConvexHullGrahamScanTBB::RunImpl() {
  const auto &data = GetInput();

  if (data.size() < 2) {
    GetOutput() = data;
    return true;
  }

  tbb::global_control gc(tbb::global_control::max_allowed_parallelism, ppc::util::GetNumThreads());

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

size_t PerepelkinIConvexHullGrahamScanTBB::FindPivotParallel(const std::vector<std::pair<double, double>> &pts) {
  auto better = [&](size_t a, size_t b) {
    if (pts[b].second < pts[a].second || (pts[b].second == pts[a].second && pts[b].first < pts[a].first)) {
      return b;
    }
    return a;
  };

  size_t pivot = tbb::parallel_reduce(tbb::blocked_range<size_t>(1, pts.size()), static_cast<size_t>(0),
                                      [&](const tbb::blocked_range<size_t> &r, size_t local_idx) {
    for (size_t i = r.begin(); i != r.end(); i++) {
      local_idx = better(local_idx, i);
    }
    return local_idx;
  }, [&](size_t a, size_t b) { return better(a, b); });

  return pivot;
}

void PerepelkinIConvexHullGrahamScanTBB::ParallelSort(std::vector<std::pair<double, double>> &data,
                                                      const std::pair<double, double> &pivot) {
  size_t n = data.size();

  if (n < 10000) {
    std::ranges::sort(data, [&](const auto &a, const auto &b) { return AngleCmp(a, b, pivot); });
    return;
  }

  tbb::parallel_sort(data.begin(), data.end(), [&](const auto &a, const auto &b) { return AngleCmp(a, b, pivot); });
}

void perepelkin_i_convex_hull_graham_scan::PerepelkinIConvexHullGrahamScanTBB::HullConstruction(
    std::vector<std::pair<double, double>> &hull, const std::vector<std::pair<double, double>> &pts,
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

double PerepelkinIConvexHullGrahamScanTBB::Orientation(const std::pair<double, double> &p,
                                                       const std::pair<double, double> &q,
                                                       const std::pair<double, double> &r) {
  double val = ((q.first - p.first) * (r.second - p.second)) - ((q.second - p.second) * (r.first - p.first));

  if (std::abs(val) < 1e-9) {
    return 0.0;
  }

  return val;
}

bool PerepelkinIConvexHullGrahamScanTBB::AngleCmp(const std::pair<double, double> &a,
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

bool PerepelkinIConvexHullGrahamScanTBB::PostProcessingImpl() {
  return true;
}

}  // namespace perepelkin_i_convex_hull_graham_scan
