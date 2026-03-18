#include "perepelkin_i_convex_hull_graham_scan/seq/include/ops_seq.hpp"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <utility>
#include <vector>

#include "perepelkin_i_convex_hull_graham_scan/common/include/common.hpp"

namespace perepelkin_i_convex_hull_graham_scan {

PerepelkinIConvexHullGrahamScanSEQ::PerepelkinIConvexHullGrahamScanSEQ(const InType &in) {
  SetTypeOfTask(GetStaticTypeOfTask());
  GetInput() = in;
  GetOutput() = std::vector<std::pair<double, double>>();
}

bool PerepelkinIConvexHullGrahamScanSEQ::ValidationImpl() {
  return GetOutput().empty();
}

bool PerepelkinIConvexHullGrahamScanSEQ::PreProcessingImpl() {
  return true;
}

bool PerepelkinIConvexHullGrahamScanSEQ::RunImpl() {
  const auto &data = GetInput();

  if (data.size() < 2) {
    GetOutput() = data;
    return true;
  }

  std::vector<std::pair<double, double>> pts = data;

  // Find pivot: lowest y, then lowest x
  auto pivot_it = std::ranges::min_element(pts, [](const auto &a, const auto &b) {
    if (a.second == b.second) {
      return a.first < b.first;
    }
    return a.second < b.second;
  });
  std::pair<double, double> pivot = *pivot_it;
  pts.erase(pivot_it);

  // Sort remaining points by polar angle around pivot (counter-clockwise)
  std::ranges::sort(pts, [&](const std::pair<double, double> &a, const std::pair<double, double> &b) {
    return AngleCmp(a, b, pivot);
  });

  // Build hull
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

double PerepelkinIConvexHullGrahamScanSEQ::Orientation(const std::pair<double, double> &p,
                                                       const std::pair<double, double> &q,
                                                       const std::pair<double, double> &r) {
  double val = ((q.first - p.first) * (r.second - p.second)) - ((q.second - p.second) * (r.first - p.first));
  if (std::abs(val) < 1e-9) {
    return 0.0;
  }
  return val;
}

bool PerepelkinIConvexHullGrahamScanSEQ::AngleCmp(const std::pair<double, double> &a,
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

bool PerepelkinIConvexHullGrahamScanSEQ::PostProcessingImpl() {
  return true;
}

}  // namespace perepelkin_i_convex_hull_graham_scan
