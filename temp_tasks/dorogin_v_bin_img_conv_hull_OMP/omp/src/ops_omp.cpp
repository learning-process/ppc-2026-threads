#include "dorogin_v_bin_img_conv_hull_OMP/omp/include/ops_omp.hpp"

#include <algorithm>
#include <array>
#include <cstddef>
#include <cstdint>
#include <ranges>
#include <utility>
#include <vector>

#include "dorogin_v_bin_img_conv_hull_OMP/common/include/common.hpp"

namespace dorogin_v_bin_img_conv_hull_omp {

namespace {

std::size_t Index(int x, int y, int w) {
  return (static_cast<std::size_t>(y) * static_cast<std::size_t>(w)) + static_cast<std::size_t>(x);
}

bool Inside(int x, int y, int w, int h) {
  return (x >= 0) && (y >= 0) && (x < w) && (y < h);
}

ComponentHull CollectOneComponent(const BinaryImage &img, std::vector<std::uint8_t> &visited, int w, int h, int start_x,
                                  int start_y, const std::array<int, 8> &dx, const std::array<int, 8> &dy) {
  ComponentHull pixels;
  pixels.reserve(64);
  std::vector<Point> stack;
  const std::size_t start_id = Index(start_x, start_y, w);
  visited[start_id] = 1U;
  stack.push_back(Point{.x = start_x, .y = start_y});

  while (!stack.empty()) {
    const Point p = stack.back();
    stack.pop_back();
    pixels.push_back(p);

    for (int dir = 0; dir < 8; ++dir) {
      const int nx = p.x + dx.at(static_cast<std::size_t>(dir));
      const int ny = p.y + dy.at(static_cast<std::size_t>(dir));
      if (!Inside(nx, ny, w, h)) {
        continue;
      }
      const std::size_t nid = Index(nx, ny, w);
      if (img.data[nid] == 0 || visited[nid] != 0U) {
        continue;
      }
      visited[nid] = 1U;
      stack.push_back(Point{.x = nx, .y = ny});
    }
  }

  return pixels;
}

std::vector<ComponentHull> CollectComponents(const BinaryImage &img) {
  const int w = img.width;
  const int h = img.height;
  const std::size_t n = static_cast<std::size_t>(w) * static_cast<std::size_t>(h);

  std::vector<std::uint8_t> visited(n, 0U);
  std::vector<ComponentHull> components;

  const std::array<int, 8> dx = {1, 1, 0, -1, -1, -1, 0, 1};
  const std::array<int, 8> dy = {0, 1, 1, 1, 0, -1, -1, -1};

  for (int yy = 0; yy < h; ++yy) {
    for (int xx = 0; xx < w; ++xx) {
      const std::size_t id = Index(xx, yy, w);
      if (img.data[id] == 0 || visited[id] != 0U) {
        continue;
      }
      ComponentHull pixels = CollectOneComponent(img, visited, w, h, xx, yy, dx, dy);
      if (!pixels.empty()) {
        components.push_back(std::move(pixels));
      }
    }
  }

  return components;
}

std::int64_t Cross(const Point &o, const Point &a, const Point &b) {
  const std::int64_t x1 = static_cast<std::int64_t>(a.x) - static_cast<std::int64_t>(o.x);
  const std::int64_t y1 = static_cast<std::int64_t>(a.y) - static_cast<std::int64_t>(o.y);
  const std::int64_t x2 = static_cast<std::int64_t>(b.x) - static_cast<std::int64_t>(o.x);
  const std::int64_t y2 = static_cast<std::int64_t>(b.y) - static_cast<std::int64_t>(o.y);
  return (x1 * y2) - (y1 * x2);
}

ComponentHull BuildHull(ComponentHull pts) {
  if (pts.size() <= 1) {
    return pts;
  }

  std::ranges::sort(pts, [](const Point &a, const Point &b) {
    if (a.x != b.x) {
      return a.x < b.x;
    }
    return a.y < b.y;
  });
  const auto uniq_result =
      std::ranges::unique(pts, [](const Point &a, const Point &b) { return (a.x == b.x) && (a.y == b.y); });
  pts.erase(uniq_result.end(), pts.end());

  if (pts.size() <= 1) {
    return pts;
  }

  ComponentHull lower;
  lower.reserve(pts.size());
  for (const auto &p : pts) {
    while (lower.size() >= 2 && Cross(lower[lower.size() - 2], lower[lower.size() - 1], p) <= 0) {
      lower.pop_back();
    }
    lower.push_back(p);
  }

  ComponentHull upper;
  upper.reserve(pts.size());
  for (std::size_t i = pts.size(); i-- > 0;) {
    const Point &p = pts[i];
    while (upper.size() >= 2 && Cross(upper[upper.size() - 2], upper[upper.size() - 1], p) <= 0) {
      upper.pop_back();
    }
    upper.push_back(p);
  }

  if (!lower.empty()) {
    lower.pop_back();
  }
  if (!upper.empty()) {
    upper.pop_back();
  }

  lower.insert(lower.end(), upper.begin(), upper.end());
  return lower;
}

}  // namespace

DoroginVBinImgConvHullOMP::DoroginVBinImgConvHullOMP(const InType &in) {
  SetTypeOfTask(GetStaticTypeOfTask());
  GetInput() = in;
  GetOutput().clear();
}

bool DoroginVBinImgConvHullOMP::ValidationImpl() {
  const auto &img = GetInput();
  if (img.width <= 0 || img.height <= 0) {
    return false;
  }
  const std::size_t required = static_cast<std::size_t>(img.width) * static_cast<std::size_t>(img.height);
  return img.data.size() == required;
}

bool DoroginVBinImgConvHullOMP::PreProcessingImpl() {
  GetOutput().clear();
  return true;
}

bool DoroginVBinImgConvHullOMP::RunImpl() {
  if (!ValidationImpl()) {
    return false;
  }

  auto components = CollectComponents(GetInput());
  OutType hulls;
  hulls.resize(components.size());

  const int count = static_cast<int>(components.size());

#pragma omp parallel for default(none) shared(components, hulls, count)
  for (int i = 0; i < count; ++i) {
    const auto idx = static_cast<std::size_t>(i);
    hulls[idx] = BuildHull(std::move(components[idx]));
  }

  GetOutput() = std::move(hulls);
  return true;
}

bool DoroginVBinImgConvHullOMP::PostProcessingImpl() {
  return true;
}

}  // namespace dorogin_v_bin_img_conv_hull_omp
