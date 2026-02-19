// tasks/peryashkin_v_binary_component_contour_processing/seq/src/ops_seq.cpp
#include "peryashkin_v_binary_component_contour_processing/seq/include/ops_seq.hpp"

#include <algorithm>
#include <array>
#include <cstddef>
#include <cstdint>
#include <queue>
#include <utility>
#include <vector>

#include "peryashkin_v_binary_component_contour_processing/common/include/common.hpp"

namespace peryashkin_v_binary_component_contour_processing {

namespace {

inline bool InBounds(int x, int y, int w, int h) {
  return x >= 0 && y >= 0 && x < w && y < h;
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

  auto uniq_it = std::ranges::unique(pts, [](const Point &a, const Point &b) { return a.x == b.x && a.y == b.y; });
  pts.erase(uniq_it.begin(), pts.end());

  if (pts.size() == 1) {
    return pts;
  }

  std::vector<Point> lower;
  lower.reserve(pts.size());
  for (const auto &pt : pts) {
    while (lower.size() >= 2 && Cross(lower[lower.size() - 2], lower[lower.size() - 1], pt) <= 0) {
      lower.pop_back();
    }
    lower.push_back(pt);
  }

  std::vector<Point> upper;
  upper.reserve(pts.size());
  for (std::size_t idx = pts.size(); idx-- > 0;) {
    const auto &pt = pts[idx];
    while (upper.size() >= 2 && Cross(upper[upper.size() - 2], upper[upper.size() - 1], pt) <= 0) {
      upper.pop_back();
    }
    upper.push_back(pt);
  }

  lower.pop_back();
  upper.pop_back();
  lower.insert(lower.end(), upper.begin(), upper.end());

  return lower;
}

inline std::size_t ToIndex(int x, int y, int w) {
  return (static_cast<std::size_t>(y) * static_cast<std::size_t>(w)) + static_cast<std::size_t>(x);
}

inline std::vector<Point> FloodFillComponent4(const BinaryImage &img, int start_x, int start_y,
                                              std::vector<std::uint8_t> &vis) {
  const int w = img.width;
  const int h = img.height;

  std::queue<Point> q;
  const std::size_t start_id = ToIndex(start_x, start_y, w);
  vis[start_id] = 1U;
  q.push(Point{.x = start_x, .y = start_y});

  std::vector<Point> pts;
  pts.reserve(128);

  static constexpr std::array<int, 4> kDx = {1, -1, 0, 0};
  static constexpr std::array<int, 4> kDy = {0, 0, 1, -1};

  while (!q.empty()) {
    const Point p = q.front();
    q.pop();
    pts.push_back(p);

    for (std::size_t dir = 0; dir < kDx.size(); ++dir) {
      const int nx = p.x + kDx[dir];
      const int ny = p.y + kDy[dir];
      if (!InBounds(nx, ny, w, h)) {
        continue;
      }
      const std::size_t nid = ToIndex(nx, ny, w);
      if (img.data[nid] == 1U && (vis[nid] == 0U)) {
        vis[nid] = 1U;
        q.push(Point{.x = nx, .y = ny});
      }
    }
  }

  return pts;
}

inline std::vector<std::vector<Point>> ExtractComponents4(const BinaryImage &img) {
  const int w = img.width;
  const int h = img.height;
  const std::size_t n = static_cast<std::size_t>(w) * static_cast<std::size_t>(h);

  std::vector<std::uint8_t> vis(n, 0U);
  std::vector<std::vector<Point>> comps;

  for (int y_pos = 0; y_pos < h; ++y_pos) {
    for (int x_pos = 0; x_pos < w; ++x_pos) {
      const std::size_t id = ToIndex(x_pos, y_pos, w);
      if (img.data[id] == 0U || (vis[id] != 0U)) {
        continue;
      }
      comps.push_back(FloodFillComponent4(img, x_pos, y_pos, vis));
    }
  }

  return comps;
}

inline OutType SolveSEQ(const BinaryImage &img) {
  auto comps = ExtractComponents4(img);
  OutType hulls;
  hulls.reserve(comps.size());
  for (auto &comp : comps) {
    hulls.push_back(ConvexHullMonotonicChain(std::move(comp)));
  }
  return hulls;
}

}  // namespace

PeryashkinVBinaryComponentContourProcessingSEQ::PeryashkinVBinaryComponentContourProcessingSEQ(const InType &in) {
  SetTypeOfTask(GetStaticTypeOfTask());
  GetInput() = in;
  GetOutput().clear();
}

bool PeryashkinVBinaryComponentContourProcessingSEQ::ValidationImpl() {
  const auto &in = GetInput();
  if (in.width <= 0 || in.height <= 0) {
    return false;
  }
  const std::size_t need = static_cast<std::size_t>(in.width) * static_cast<std::size_t>(in.height);
  return in.data.size() == need;
}

bool PeryashkinVBinaryComponentContourProcessingSEQ::PreProcessingImpl() {
  local_out_.clear();
  return true;
}

bool PeryashkinVBinaryComponentContourProcessingSEQ::RunImpl() {
  if (!ValidationImpl()) {
    return false;
  }

  local_out_ = SolveSEQ(GetInput());
  return true;
}

bool PeryashkinVBinaryComponentContourProcessingSEQ::PostProcessingImpl() {
  GetOutput() = local_out_;
  return true;
}

}  // namespace peryashkin_v_binary_component_contour_processing
