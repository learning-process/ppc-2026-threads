#include "yurkin_g_graham_scan/seq/include/ops_seq.hpp"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <vector>
#include <ranges>

#include "yurkin_g_graham_scan/common/include/common.hpp"

namespace yurkin_g_graham_scan {

YurkinGGrahamScanSEQ::YurkinGGrahamScanSEQ(const InType &in) {
  SetTypeOfTask(GetStaticTypeOfTask());
  GetInput() = in;
  GetOutput().clear();
}

namespace {  // anonymous namespace for internal helpers
long double Cross(const Point &o, const Point &a, const Point &b) {
  return (static_cast<long double>(a.x - o.x) * static_cast<long double>(b.y - o.y))
         - (static_cast<long double>(a.y - o.y) * static_cast<long double>(b.x - o.x));
}
}  // namespace

bool YurkinGGrahamScanSEQ::ValidationImpl() {
  const auto &pts = GetInput();
  return !pts.empty();
}

bool YurkinGGrahamScanSEQ::PreProcessingImpl() {
  auto &pts = GetInput();
  if (pts.empty()) return true;

  // Sort using ranges
  std::ranges::sort(pts, [](const Point &a, const Point &b) {
    if (a.x != b.x) return a.x < b.x;
    return a.y < b.y;
  });

  // Remove exact duplicates explicitly (stable, clear, avoids std::unique warning)
  std::vector<Point> tmp;
  tmp.reserve(pts.size());
  for (const auto &p : pts) {
    if (tmp.empty() || tmp.back().x != p.x || tmp.back().y != p.y) {
      tmp.push_back(p);
    }
  }
  pts.swap(tmp);

  return !pts.empty();
}

bool YurkinGGrahamScanSEQ::RunImpl() {
  const InType &pts_in = GetInput();
  const std::size_t n = pts_in.size();
  if (n == 0) {
    GetOutput().clear();
    return true;
  }
  if (n == 1) {
    GetOutput() = pts_in;
    return true;
  }

  InType pts = pts_in;
  std::ranges::sort(pts, [](const Point &a, const Point &b) {
    if (a.x != b.x) return a.x < b.x;
    return a.y < b.y;
  });

  OutType lower;
  lower.reserve(pts.size());
  for (const auto &p : pts) {
    while (lower.size() >= 2 &&
           Cross(lower[lower.size() - 2], lower[lower.size() - 1], p) <= 0) {
      lower.pop_back();
    }
    lower.push_back(p);
  }

  OutType upper;
  upper.reserve(pts.size());
  for (const auto &p : std::ranges::reverse_view(pts)) {
    while (upper.size() >= 2 &&
           Cross(upper[upper.size() - 2], upper[upper.size() - 1], p) <= 0) {
      upper.pop_back();
    }
    upper.push_back(p);
  }

  OutType hull;
  hull.reserve(lower.size() + upper.size());
  for (const auto &pt : lower) hull.push_back(pt);
  for (std::size_t i = 1; i + 1 < upper.size(); ++i) hull.push_back(upper[i]);

  GetOutput() = hull;
  return true;
}

bool YurkinGGrahamScanSEQ::PostProcessingImpl() {
  if (GetInput().empty()) return true;
  return !GetOutput().empty();
}

}  // namespace yurkin_g_graham_scan
