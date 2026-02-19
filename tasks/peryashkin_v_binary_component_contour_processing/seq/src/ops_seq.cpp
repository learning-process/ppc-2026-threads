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

  std::sort(pts.begin(), pts.end(),
            [](const Point &a, const Point &b) { return (a.x < b.x) || ((a.x == b.x) && (a.y < b.y)); });

  pts.erase(
      std::unique(pts.begin(), pts.end(), [](const Point &a, const Point &b) { return (a.x == b.x) && (a.y == b.y); }),
      pts.end());

  if (pts.size() == 1) {
    return pts;
  }

  std::vector<Point> lower;
  lower.reserve(pts.size());
  for (const auto &p : pts) {
    while (lower.size() >= 2 && Cross(lower[lower.size() - 2], lower[lower.size() - 1], p) <= 0) {
      lower.pop_back();
    }
    lower.push_back(p);
  }

  std::vector<Point> upper;
  upper.reserve(pts.size());
  for (std::size_t i = pts.size(); i-- > 0;) {
    const auto &p = pts[i];
    while (upper.size() >= 2 && Cross(upper[upper.size() - 2], upper[upper.size() - 1], p) <= 0) {
      upper.pop_back();
    }
    upper.push_back(p);
  }

  lower.pop_back();
  upper.pop_back();
  lower.insert(lower.end(), upper.begin(), upper.end());

  return lower;
}

inline std::vector<std::vector<Point>> ExtractComponents4(const BinaryImage &img) {
  const int w = img.width;
  const int h = img.height;
  const std::size_t n = static_cast<std::size_t>(w) * static_cast<std::size_t>(h);

  std::vector<std::uint8_t> vis(n, 0);
  std::vector<std::vector<Point>> comps;

  auto idx = [w](int x, int y) -> std::size_t {
    return (static_cast<std::size_t>(y) * static_cast<std::size_t>(w)) + static_cast<std::size_t>(x);
  };

  static constexpr std::array<int, 4> kDx = {1, -1, 0, 0};
  static constexpr std::array<int, 4> kDy = {0, 0, 1, -1};

  std::queue<Point> q;
  for (int yy = 0; yy < h; ++yy) {
    for (int xx = 0; xx < w; ++xx) {
      const std::size_t id = idx(xx, yy);
      if ((img.data[id] == 0) || (vis[id] != 0u)) {
        continue;
      }

      vis[id] = 1;
      q.push(Point{xx, yy});

      std::vector<Point> pts;
      pts.reserve(128);

      while (!q.empty()) {
        const Point p = q.front();
        q.pop();
        pts.push_back(p);

        for (int dir = 0; dir < 4; ++dir) {
          const int nx = p.x + kDx.at(static_cast<std::size_t>(dir));
          const int ny = p.y + kDy.at(static_cast<std::size_t>(dir));
          if (!InBounds(nx, ny, w, h)) {
            continue;
          }

          const std::size_t nid = idx(nx, ny);
          if ((img.data[nid] == 1) && (vis[nid] == 0u)) {
            vis[nid] = 1;
            q.push(Point{nx, ny});
          }
        }
      }

      comps.push_back(std::move(pts));
    }
  }

  return comps;
}

inline OutType SolveSEQ(const BinaryImage &img) {
  auto comps = ExtractComponents4(img);
  OutType hulls;
  hulls.reserve(comps.size());
  for (auto &c : comps) {
    hulls.push_back(ConvexHullMonotonicChain(std::move(c)));
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
