#include "dorogin_v_bin_img_conv_hull/seq/include/ops_seq.hpp"

#include <algorithm>
#include <array>
#include <cstddef>
#include <cstdint>
#include <ranges>
#include <stack>
#include <utility>
#include <vector>
#include <functional>

#include "dorogin_v_bin_img_conv_hull/common/include/common.hpp"

namespace dorogin_v_bin_img_conv_hull {

namespace {

constexpr uint8_t kThreshold = 128;

inline bool IsInside(int col, int row, int width, int height) {
  return col >= 0 && row >= 0 && col < width && row < height;
}

int64_t Cross(const Point &a, const Point &b, const Point &c) {
  int64_t x1 = b.x - a.x;
  int64_t y1 = b.y - a.y;
  int64_t x2 = c.x - a.x;
  int64_t y2 = c.y - a.y;
  return (x1 * y2) - (y1 * x2);
}

}  // namespace

DoroginVBinImgConvHullSeq::DoroginVBinImgConvHullSeq(const InType &in) : w_(in) {
  SetTypeOfTask(GetStaticTypeOfTask());
  GetInput() = in;
}

bool DoroginVBinImgConvHullSeq::ValidationImpl() {
  const auto &in = GetInput();
  if (in.width <= 0 || in.height <= 0) {
    return false;
  }
  return in.pixels.size() == static_cast<size_t>(in.width) * static_cast<size_t>(in.height);
}

bool DoroginVBinImgConvHullSeq::PreProcessingImpl() {
  w_ = GetInput();
  ThresholdImage();
  return true;
}

bool DoroginVBinImgConvHullSeq::RunImpl() {
  FindComponents();

  w_.convex_hulls.clear();

  for (const auto &comp : w_.components) {
    if (comp.empty()) {
      continue;
    }

    if (comp.size() < 3) {
      w_.convex_hulls.push_back(comp);
    } else {
      w_.convex_hulls.push_back(BuildHull(comp));
    }
  }

  GetOutput() = w_;
  return true;
}

bool DoroginVBinImgConvHullSeq::PostProcessingImpl() {
  return true;
}

size_t DoroginVBinImgConvHullSeq::Index(int col, int row, int width) {
  return (static_cast<size_t>(row) * static_cast<size_t>(width)) + static_cast<size_t>(col);
}

void DoroginVBinImgConvHullSeq::ThresholdImage() {
  std::ranges::transform(w_.pixels, w_.pixels.begin(), [](uint8_t p) { return p > kThreshold ? 255 : 0; });
}

void DoroginVBinImgConvHullSeq::ExploreComponent(int start_col, int start_row, int width, int height,
                                                 std::vector<bool> &visited, std::vector<Point> &component) {
  std::stack<Point> stack;
  stack.emplace(start_col, start_row);

  visited[Index(start_col, start_row, width)] = true;

  const std::array<int, 4> dx{1, -1, 0, 0};
  const std::array<int, 4> dy{0, 0, 1, -1};

  while (!stack.empty()) {
    Point p = stack.top();
    stack.pop();

    component.push_back(p);

    for (size_t dir = 0; dir < dx.size(); ++dir) {
      const int nx = p.x + dx.at(dir);
      const int ny = p.y + dy.at(dir);

      if (!IsInside(nx, ny, width, height)) {
        continue;
      }

      size_t idx = Index(nx, ny, width);
      if (visited[idx]) {
        continue;
      }
      if (w_.pixels[idx] == 0) {
        continue;
      }

      visited[idx] = true;
      stack.emplace(nx, ny);
    }
  }
}

void DoroginVBinImgConvHullSeq::FindComponents() {
  const int width = w_.width;
  const int height = w_.height;

  std::vector<bool> visited(static_cast<size_t>(width) * height, false);
  w_.components.clear();

  for (int row = 0; row < height; ++row) {
    for (int col = 0; col < width; ++col) {
      size_t idx = Index(col, row, width);

      if (w_.pixels[idx] == 0 || visited[idx]) {
        continue;
      }

      std::vector<Point> comp;
      ExploreComponent(col, row, width, height, visited, comp);

      w_.components.push_back(std::move(comp));
    }
  }
}

std::vector<Point> DoroginVBinImgConvHullSeq::BuildHull(const std::vector<Point> &points) {
  std::vector<Point> pts = points;

  std::ranges::sort(pts, std::less<>{});
  auto unique_end = std::ranges::unique(pts).begin();
  pts.erase(unique_end, pts.end());

  if (pts.size() < 3) {
    return pts;
  }

  std::vector<Point> hull;

  for (const auto &p : pts) {
    while (hull.size() >= 2 && Cross(hull[hull.size() - 2], hull.back(), p) <= 0) {
      hull.pop_back();
    }
    hull.push_back(p);
  }

  size_t lsize = hull.size();
  for (int i = static_cast<int>(pts.size()) - 2; i >= 0; --i) {
    while (hull.size() > lsize && Cross(hull[hull.size() - 2], hull.back(), pts[i]) <= 0) {
      hull.pop_back();
    }
    hull.push_back(pts[i]);
  }

  hull.pop_back();
  return hull;
}

}  // namespace dorogin_v_bin_img_conv_hull
