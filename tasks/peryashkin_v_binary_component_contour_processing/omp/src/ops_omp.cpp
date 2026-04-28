#include "peryashkin_v_binary_component_contour_processing/omp/include/ops_omp.hpp"

#include <algorithm>
#include <cstddef>
#include <cstdint>
#include <queue>
#include <utility>
#include <vector>

#include "peryashkin_v_binary_component_contour_processing/common/include/common.hpp"

namespace peryashkin_v_binary_component_contour_processing {

namespace {

inline bool InBounds(int x, int y, int w, int h) {
  return (x >= 0) && (y >= 0) && (x < w) && (y < h);
}

inline std::int64_t Cross(const Point &o, const Point &a, const Point &b) {
  const std::int64_t x1 = static_cast<std::int64_t>(a.x) - static_cast<std::int64_t>(o.x);
  const std::int64_t y1 = static_cast<std::int64_t>(a.y) - static_cast<std::int64_t>(o.y);
  const std::int64_t x2 = static_cast<std::int64_t>(b.x) - static_cast<std::int64_t>(o.x);
  const std::int64_t y2 = static_cast<std::int64_t>(b.y) - static_cast<std::int64_t>(o.y);
  return (x1 * y2) - (y1 * x2);
}

inline std::vector<Point> ConvexHullMonotonicChain(std::vector<Point> pts) {
  if (pts.empty()) {
    return {};
  }

  std::ranges::sort(pts, [](const Point &a, const Point &b) { return (a.x < b.x) || ((a.x == b.x) && (a.y < b.y)); });

  const auto new_end =
      std::ranges::unique(pts, [](const Point &a, const Point &b) { return (a.x == b.x) && (a.y == b.y); }).begin();
  pts.erase(new_end, pts.end());

  if (pts.size() == 1) {
    return pts;
  }

  std::vector<Point> lower;
  lower.reserve(pts.size());
  for (const auto &p : pts) {
    while ((lower.size() >= 2) && (Cross(lower[lower.size() - 2], lower[lower.size() - 1], p) <= 0)) {
      lower.pop_back();
    }
    lower.push_back(p);
  }

  std::vector<Point> upper;
  upper.reserve(pts.size());
  for (std::size_t i = pts.size(); i-- > 0;) {
    const auto &p = pts[i];
    while ((upper.size() >= 2) && (Cross(upper[upper.size() - 2], upper[upper.size() - 1], p) <= 0)) {
      upper.pop_back();
    }
    upper.push_back(p);
  }

  lower.pop_back();
  upper.pop_back();
  lower.insert(lower.end(), upper.begin(), upper.end());
  return lower;
}

inline std::size_t Index2D(int x, int y, int w) {
  return (static_cast<std::size_t>(y) * static_cast<std::size_t>(w)) + static_cast<std::size_t>(x);
}

inline void TryPushNeighbor(const BinaryImage &img, std::vector<std::uint8_t> &vis, std::queue<Point> &q, int w, int h,
                            int nx, int ny) {
  if (!InBounds(nx, ny, w, h)) {
    return;
  }
  const std::size_t nid = Index2D(nx, ny, w);
  if ((img.data[nid] == 1) && (vis[nid] == 0U)) {
    vis[nid] = 1U;
    q.push(Point{.x = nx, .y = ny});
  }
}

inline std::vector<Point> BfsComponent4(const BinaryImage &img, std::vector<std::uint8_t> &vis, int w, int h, int sx,
                                        int sy) {
  std::queue<Point> q;
  const std::size_t start_id = Index2D(sx, sy, w);
  vis[start_id] = 1U;
  q.push(Point{.x = sx, .y = sy});

  std::vector<Point> pts;
  pts.reserve(128);

  while (!q.empty()) {
    const Point p = q.front();
    q.pop();
    pts.push_back(p);

    TryPushNeighbor(img, vis, q, w, h, p.x + 1, p.y);
    TryPushNeighbor(img, vis, q, w, h, p.x - 1, p.y);
    TryPushNeighbor(img, vis, q, w, h, p.x, p.y + 1);
    TryPushNeighbor(img, vis, q, w, h, p.x, p.y - 1);
  }

  return pts;
}

inline std::vector<std::vector<Point>> ExtractComponents4(const BinaryImage &img) {
  const int w = img.width;
  const int h = img.height;
  const std::size_t n = static_cast<std::size_t>(w) * static_cast<std::size_t>(h);

  std::vector<std::uint8_t> vis(n, 0U);
  std::vector<std::vector<Point>> comps;

  for (int yy = 0; yy < h; ++yy) {
    for (int xx = 0; xx < w; ++xx) {
      const std::size_t id = Index2D(xx, yy, w);
      if ((img.data[id] == 0) || (vis[id] != 0U)) {
        continue;
      }
      comps.push_back(BfsComponent4(img, vis, w, h, xx, yy));
    }
  }

  return comps;
}

inline OutType SolveOMP(const BinaryImage &img) {
  auto comps = ExtractComponents4(img);
  OutType hulls;
  hulls.resize(comps.size());

  const int n = static_cast<int>(comps.size());

#pragma omp parallel for default(none) shared(comps, hulls, n) schedule(static)
  for (int i = 0; i < n; ++i) {
    const auto idx = static_cast<std::size_t>(i);
    hulls[idx] = ConvexHullMonotonicChain(std::move(comps[idx]));
  }

  return hulls;
}

}  // namespace

PeryashkinVBinaryComponentContourProcessingOMP::PeryashkinVBinaryComponentContourProcessingOMP(const InType &in) {
  SetTypeOfTask(GetStaticTypeOfTask());
  GetInput() = in;
  GetOutput().clear();
}

bool PeryashkinVBinaryComponentContourProcessingOMP::ValidationImpl() {
  const auto &in = GetInput();
  if ((in.width <= 0) || (in.height <= 0)) {
    return false;
  }
  const std::size_t need = static_cast<std::size_t>(in.width) * static_cast<std::size_t>(in.height);
  return in.data.size() == need;
}

bool PeryashkinVBinaryComponentContourProcessingOMP::PreProcessingImpl() {
  local_out_.clear();
  return true;
}

bool PeryashkinVBinaryComponentContourProcessingOMP::RunImpl() {
  if (!ValidationImpl()) {
    return false;
  }
  local_out_ = SolveOMP(GetInput());
  return true;
}

bool PeryashkinVBinaryComponentContourProcessingOMP::PostProcessingImpl() {
  GetOutput() = local_out_;
  return true;
}

}  // namespace peryashkin_v_binary_component_contour_processing
