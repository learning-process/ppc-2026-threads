#include "shkenev_i_constra_hull_for_binary_image/tbb/include/ops_tbb.hpp"

#include <tbb/blocked_range.h>
#include <tbb/parallel_for.h>

#include <algorithm>
#include <array>
#include <cstddef>
#include <cstdint>
#include <queue>
#include <ranges>
#include <utility>
#include <vector>

#include "shkenev_i_constra_hull_for_binary_image/common/include/common.hpp"

namespace shkenev_i_constra_hull_for_binary_image {

namespace {

// 12
constexpr uint8_t kThreshold = 128;

constexpr std::array<std::pair<int, int>, 4> kDirs = {{{1, 0}, {-1, 0}, {0, 1}, {0, -1}}};

int64_t Cross(const Point &a, const Point &b, const Point &c) {
  return (static_cast<int64_t>(b.x - a.x) * (c.y - b.y)) - (static_cast<int64_t>(b.y - a.y) * (c.x - b.x));
}

inline bool InBounds(int x, int y, int w, int h) {
  return x >= 0 && x < w && y >= 0 && y < h;
}

}  // namespace

ShkenevIConstrHullTBB::ShkenevIConstrHullTBB(const InType &in) : work_(in) {
  SetTypeOfTask(GetStaticTypeOfTask());
  GetInput() = in;
}

bool ShkenevIConstrHullTBB::ValidationImpl() {
  const auto &in = GetInput();
  return in.width > 0 && in.height > 0 &&
         in.pixels.size() == static_cast<size_t>(in.width) * static_cast<size_t>(in.height);
}

bool ShkenevIConstrHullTBB::PreProcessingImpl() {
  work_ = GetInput();
  work_.components.clear();
  work_.convex_hulls.clear();

  ThresholdImage();
  return true;
}

void ShkenevIConstrHullTBB::ThresholdImage() {
  auto &p = work_.pixels;

  tbb::parallel_for(tbb::blocked_range<size_t>(0, p.size()), [&](const auto &r) {
    for (size_t i = r.begin(); i < r.end(); ++i) {
      p[i] = (p[i] > kThreshold) ? 255 : 0;
    }
  });
}

void ShkenevIConstrHullTBB::AddNeighbors(const Point &current, int width, int height, std::vector<uint8_t> &visited,
                                         std::queue<Point> &queue) {
  for (auto [dx, dy] : kDirs) {
    int neighbor_x = current.x + dx;
    int neighbor_y = current.y + dy;

    if (!InBounds(neighbor_x, neighbor_y, width, height)) {
      continue;
    }

    size_t neighbor_idx = Index(neighbor_x, neighbor_y, width);

    if (visited[neighbor_idx] != 0U || work_.pixels[neighbor_idx] == 0) {
      continue;
    }

    visited[neighbor_idx] = 1;
    queue.emplace(neighbor_x, neighbor_y);
  }
}

std::vector<Point> ShkenevIConstrHullTBB::ExtractComponent(int start_x, int start_y, int width, int height,
                                                           std::vector<uint8_t> &visited) {
  std::vector<Point> component;
  std::queue<Point> queue;
  size_t start_idx = Index(start_x, start_y, width);

  queue.emplace(start_x, start_y);
  visited[start_idx] = 1;

  while (!queue.empty()) {
    auto current = queue.front();
    queue.pop();
    component.push_back(current);
    AddNeighbors(current, width, height, visited, queue);
  }

  return component;
}

void ShkenevIConstrHullTBB::FindComponents() {
  const int width = work_.width;
  const int height = work_.height;
  std::vector<uint8_t> visited(static_cast<size_t>(width) * height, 0);
  work_.components.clear();

  for (int y_coord = 0; y_coord < height; ++y_coord) {
    for (int x_coord = 0; x_coord < width; ++x_coord) {
      size_t idx = Index(x_coord, y_coord, width);

      if (visited[idx] != 0U || work_.pixels[idx] == 0) {
        continue;
      }

      work_.components.push_back(ExtractComponent(x_coord, y_coord, width, height, visited));
    }
  }
}

bool ShkenevIConstrHullTBB::RunImpl() {
  FindComponents();

  auto &comps = work_.components;
  auto &hulls = work_.convex_hulls;

  hulls.resize(comps.size());

  tbb::parallel_for(tbb::blocked_range<size_t>(0, comps.size()), [&](const auto &r) {
    for (size_t i = r.begin(); i < r.end(); ++i) {
      const auto &comp = comps[i];

      if (comp.size() <= 2) {
        hulls[i] = comp;
      } else {
        hulls[i] = BuildHull(comp);
      }
    }
  });

  GetOutput() = work_;
  return true;
}

std::vector<Point> ShkenevIConstrHullTBB::BuildHull(const std::vector<Point> &points) {
  if (points.size() <= 2) {
    return points;
  }

  std::vector<Point> pts = points;

  std::ranges::sort(pts, [](auto &a, auto &b) { return (a.x != b.x) ? (a.x < b.x) : (a.y < b.y); });

  auto [first, last] = std::ranges::unique(pts);
  pts.erase(first, last);

  if (pts.size() <= 2) {
    return pts;
  }

  std::vector<Point> lower;
  std::vector<Point> upper;
  lower.reserve(pts.size());
  upper.reserve(pts.size());

  for (auto &point : pts) {
    while (lower.size() >= 2 && Cross(lower[lower.size() - 2], lower.back(), point) <= 0) {
      lower.pop_back();
    }
    lower.push_back(point);
  }

  for (auto &point : std::ranges::reverse_view(pts)) {
    while (upper.size() >= 2 && Cross(upper[upper.size() - 2], upper.back(), point) <= 0) {
      upper.pop_back();
    }
    upper.push_back(point);
  }

  lower.pop_back();
  upper.pop_back();

  lower.insert(lower.end(), upper.begin(), upper.end());
  return lower;
}

size_t ShkenevIConstrHullTBB::Index(int x, int y, int width) {
  return (static_cast<size_t>(y) * static_cast<size_t>(width)) + static_cast<size_t>(x);
}

bool ShkenevIConstrHullTBB::PostProcessingImpl() {
  return true;
}

}  // namespace shkenev_i_constra_hull_for_binary_image
