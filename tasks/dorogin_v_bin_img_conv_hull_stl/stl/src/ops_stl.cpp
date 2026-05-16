#include "dorogin_v_bin_img_conv_hull_stl/stl/include/ops_stl.hpp"

#include <algorithm>
#include <cstddef>
#include <cstdint>
#include <queue>
#include <ranges>
#include <utility>
#include <vector>

#include "dorogin_v_bin_img_conv_hull_stl/common/include/common.hpp"

namespace nesterov_a_test_task_threads {
namespace {
inline bool InBounds(const int x, const int y, const int w, const int h) {
  return (x >= 0) && (y >= 0) && (x < w) && (y < h);
}

inline std::size_t Idx(const int x, const int y, const int w) {
  return (static_cast<std::size_t>(y) * static_cast<std::size_t>(w)) + static_cast<std::size_t>(x);
}

inline std::int64_t Cross(const Point &o, const Point &a, const Point &b) {
  return (static_cast<std::int64_t>(a.x - o.x) * static_cast<std::int64_t>(b.y - o.y)) -
         (static_cast<std::int64_t>(a.y - o.y) * static_cast<std::int64_t>(b.x - o.x));
}

inline std::vector<Point> ConvexHullMonotonicChain(std::vector<Point> pts) {
  if (pts.empty()) {
    return {};
  }
  std::ranges::sort(pts, [](const Point &a, const Point &b) { return (a.x < b.x) || ((a.x == b.x) && (a.y < b.y)); });
  const auto uniq =
      std::ranges::unique(pts, [](const Point &a, const Point &b) { return (a.x == b.x) && (a.y == b.y); });
  pts.erase(uniq.begin(), pts.end());
  if (pts.size() <= 1) {
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

inline void TryPush4(const BinaryImage &img, const int w, const int h, const int nx, const int ny,
                     std::vector<std::uint8_t> &vis, std::queue<Point> &q) {
  if (!InBounds(nx, ny, w, h)) {
    return;
  }
  const std::size_t nid = Idx(nx, ny, w);
  if ((img.data[nid] == 1U) && (vis[nid] == 0U)) {
    vis[nid] = 1U;
    q.push(Point{.x = nx, .y = ny});
  }
}

inline std::vector<Point> BfsComponent4(const BinaryImage &img, const int w, const int h, const int sx, const int sy,
                                        std::vector<std::uint8_t> &vis) {
  std::vector<Point> pts;
  std::queue<Point> q;
  const std::size_t sid = Idx(sx, sy, w);
  vis[sid] = 1U;
  q.push(Point{.x = sx, .y = sy});
  while (!q.empty()) {
    const Point p = q.front();
    q.pop();
    pts.push_back(p);
    TryPush4(img, w, h, p.x + 1, p.y, vis, q);
    TryPush4(img, w, h, p.x - 1, p.y, vis, q);
    TryPush4(img, w, h, p.x, p.y + 1, vis, q);
    TryPush4(img, w, h, p.x, p.y - 1, vis, q);
  }
  return pts;
}

inline std::vector<std::vector<Point>> ExtractComponents4(const BinaryImage &img) {
  std::vector<std::vector<Point>> comps;
  if ((img.width <= 0) || (img.height <= 0)) {
    return comps;
  }
  const std::size_t n = static_cast<std::size_t>(img.width) * static_cast<std::size_t>(img.height);
  std::vector<std::uint8_t> vis(n, 0U);
  for (int row = 0; row < img.height; ++row) {
    for (int col = 0; col < img.width; ++col) {
      const std::size_t id = Idx(col, row, img.width);
      if ((img.data[id] == 0U) || (vis[id] != 0U)) {
        continue;
      }
      comps.push_back(BfsComponent4(img, img.width, img.height, col, row, vis));
    }
  }
  return comps;
}

inline bool ValidateInput(const BinaryImage &img) {
  if ((img.width <= 0) || (img.height <= 0)) {
    return false;
  }
  const std::size_t need = static_cast<std::size_t>(img.width) * static_cast<std::size_t>(img.height);
  return img.data.size() == need;
}

inline OutType SolveSTL(const BinaryImage &img) {
  auto comps = ExtractComponents4(img);
  OutType hulls;
  hulls.resize(comps.size());
  std::vector<std::size_t> indices(comps.size());
  for (std::size_t i = 0; i < indices.size(); ++i) {
    indices[i] = i;
  }
  for (const std::size_t i : indices) {
    hulls[i] = ConvexHullMonotonicChain(std::move(comps[i]));
  }
  return hulls;
}

}  // namespace

NesterovATestTaskSTL::NesterovATestTaskSTL(const InType &in) {
  SetTypeOfTask(GetStaticTypeOfTask());
  GetInput() = in;
  GetOutput().clear();
}

bool NesterovATestTaskSTL::ValidationImpl() {
  return ValidateInput(GetInput());
}

bool NesterovATestTaskSTL::PreProcessingImpl() {
  local_out_.clear();
  return true;
}

bool NesterovATestTaskSTL::RunImpl() {
  if (!ValidationImpl()) {
    return false;
  }

  local_out_ = SolveSTL(GetInput());
  return true;
}

bool NesterovATestTaskSTL::PostProcessingImpl() {
  GetOutput() = local_out_;
  return true;
}

}  // namespace nesterov_a_test_task_threads
