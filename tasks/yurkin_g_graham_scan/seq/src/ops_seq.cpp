#include "yurkin_g_graham_scan/seq/include/ops_seq.hpp"

#include <algorithm>
#include <cmath>
#include <set>
#include <vector>

#include "util/include/util.hpp"
#include "yurkin_g_graham_scan/common/include/common.hpp"

namespace yurkin_g_graham_scan {

YurkinGGrahamScanSEQ::YurkinGGrahamScanSEQ(const InType &in) {
  SetTypeOfTask(GetStaticTypeOfTask());
  GetInput() = in;
  GetOutput().clear();
}

static long double cross(const Point &O, const Point &A, const Point &B) {
  return (long double)(A.x - O.x) * (long double)(B.y - O.y) - (long double)(A.y - O.y) * (long double)(B.x - O.x);
}

bool YurkinGGrahamScanSEQ::ValidationImpl() {
  const auto &pts = GetInput();
  return !pts.empty();
}

bool YurkinGGrahamScanSEQ::PreProcessingImpl() {
  auto &pts = GetInput();
  if (pts.empty()) {
    return true;
  }
  std::sort(pts.begin(), pts.end(), [](const Point &a, const Point &b) {
    if (a.x != b.x) {
      return a.x < b.x;
    }
    return a.y < b.y;
  });
  pts.erase(
      std::unique(pts.begin(), pts.end(), [](const Point &a, const Point &b) { return a.x == b.x && a.y == b.y; }),
      pts.end());
  return !pts.empty();
}

bool YurkinGGrahamScanSEQ::RunImpl() {
  const InType &pts_in = GetInput();
  const size_t n = pts_in.size();
  if (n == 0) {
    GetOutput().clear();
    return true;
  }
  if (n == 1) {
    GetOutput() = pts_in;
    return true;
  }

  InType pts = pts_in;
  std::sort(pts.begin(), pts.end(), [](const Point &a, const Point &b) {
    if (a.x != b.x) {
      return a.x < b.x;
    }
    return a.y < b.y;
  });

  OutType lower;
  lower.reserve(pts.size());
  for (const auto &p : pts) {
    while (lower.size() >= 2 && cross(lower[lower.size() - 2], lower[lower.size() - 1], p) <= 0) {
      lower.pop_back();
    }
    lower.push_back(p);
  }

  OutType upper;
  upper.reserve(pts.size());
  for (auto it = pts.rbegin(); it != pts.rend(); ++it) {
    const auto &p = *it;
    while (upper.size() >= 2 && cross(upper[upper.size() - 2], upper[upper.size() - 1], p) <= 0) {
      upper.pop_back();
    }
    upper.push_back(p);
  }

  OutType hull;
  hull.reserve(lower.size() + upper.size());
  for (size_t i = 0; i < lower.size(); ++i) {
    hull.push_back(lower[i]);
  }
  for (size_t i = 1; i + 1 < upper.size(); ++i) {
    hull.push_back(upper[i]);
  }

  GetOutput() = hull;
  return true;
}

bool YurkinGGrahamScanSEQ::PostProcessingImpl() {
  if (GetInput().empty()) {
    return true;
  }
  return !GetOutput().empty();
}

}  // namespace yurkin_g_graham_scan
