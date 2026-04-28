#include "shkenev_i_constra_hull_for_binary_image/omp/include/ops_omp.hpp"

#include <omp.h>

#include <algorithm>
#include <array>
#include <cstddef>
#include <cstdint>
#include <queue>
#include <ranges>
#include <utility>
#include <vector>

#include "shkenev_i_constra_hull_for_binary_image/common/include/common.hpp"
#include "util/include/util.hpp"

namespace shkenev_i_constra_hull_for_binary_image {

namespace {

constexpr uint8_t kThreshold = 128;
constexpr std::array<std::pair<int, int>, 4> kDirs = {std::make_pair(1, 0), std::make_pair(-1, 0), std::make_pair(0, 1),
                                                      std::make_pair(0, -1)};

int64_t Cross(const Point &a, const Point &b, const Point &c) {
  int64_t abx = static_cast<int64_t>(b.x) - static_cast<int64_t>(a.x);
  int64_t aby = static_cast<int64_t>(b.y) - static_cast<int64_t>(a.y);
  int64_t bcx = static_cast<int64_t>(c.x) - static_cast<int64_t>(b.x);
  int64_t bcy = static_cast<int64_t>(c.y) - static_cast<int64_t>(b.y);
  return (abx * bcy) - (aby * bcx);
}

bool IsForeground(uint8_t pixel) {
  return pixel > kThreshold;
}

bool IsInBounds(int x, int y, int width, int height) {
  return x >= 0 && x < width && y >= 0 && y < height;
}

}  // namespace

ShkenevIConstrHullOMP::ShkenevIConstrHullOMP(const InType &in) : work_(in) {
  SetTypeOfTask(GetStaticTypeOfTask());
  GetInput() = in;
}

bool ShkenevIConstrHullOMP::ValidationImpl() {
  const auto &in = GetInput();
  bool dims = in.width > 0 && in.height > 0;
  bool size = in.pixels.size() == static_cast<size_t>(in.width) * static_cast<size_t>(in.height);
  return dims && size;
}

bool ShkenevIConstrHullOMP::PreProcessingImpl() {
  work_ = GetInput();
  ThresholdImage();
  return true;
}

bool ShkenevIConstrHullOMP::RunImpl() {
  FindComponents();

  work_.convex_hulls.clear();
  work_.convex_hulls.resize(work_.components.size());

  auto &components = work_.components;
  auto &convex_hulls = work_.convex_hulls;

#pragma omp parallel for default(none) shared(components, convex_hulls) num_threads(ppc::util::GetNumThreads())
  for (std::size_t i = 0; i < components.size(); ++i) {
    const auto &comp = components[i];
    if (comp.empty()) {
      continue;
    }
    if (comp.size() <= 2) {
      convex_hulls[i] = comp;
    } else {
      convex_hulls[i] = BuildHull(comp);
    }
  }

  GetOutput() = work_;
  return true;
}

bool ShkenevIConstrHullOMP::PostProcessingImpl() {
  return true;
}

size_t ShkenevIConstrHullOMP::Index(int x, int y, int width) {
  return (static_cast<size_t>(y) * static_cast<size_t>(width)) + static_cast<size_t>(x);
}

void ShkenevIConstrHullOMP::ThresholdImage() {
  auto &pixels = work_.pixels;

  for (auto &pixel : pixels) {
    pixel = IsForeground(pixel) ? static_cast<uint8_t>(255) : static_cast<uint8_t>(0);
  }
}

void ShkenevIConstrHullOMP::ExploreComponent(int start_col, int start_row, int width, int height,
                                             std::vector<uint8_t> &visited, std::vector<Point> &component) {
  std::queue<Point> queue;
  queue.emplace(start_col, start_row);
  visited[Index(start_col, start_row, width)] = 1;

  while (!queue.empty()) {
    Point current = queue.front();
    queue.pop();
    component.push_back(current);

    for (auto [dx, dy] : kDirs) {
      int next_x = current.x + dx;
      int next_y = current.y + dy;

      if (!IsInBounds(next_x, next_y, width, height)) {
        continue;
      }

      size_t next_idx = Index(next_x, next_y, width);
      if (visited[next_idx] != 0 || work_.pixels[next_idx] == 0) {
        continue;
      }

      visited[next_idx] = 1;
      queue.emplace(next_x, next_y);
    }
  }
}

void ShkenevIConstrHullOMP::FindComponents() {
  int width = work_.width;
  int height = work_.height;

  std::vector<uint8_t> visited(static_cast<size_t>(width) * static_cast<size_t>(height), 0);
  work_.components.clear();

  for (int row = 0; row < height; ++row) {
    for (int col = 0; col < width; ++col) {
      size_t idx = Index(col, row, width);
      if (work_.pixels[idx] == 0 || visited[idx] != 0) {
        continue;
      }

      std::vector<Point> component;
      ExploreComponent(col, row, width, height, visited, component);

      if (!component.empty()) {
        work_.components.push_back(std::move(component));
      }
    }
  }
}

std::vector<Point> ShkenevIConstrHullOMP::BuildHull(const std::vector<Point> &points) {
  if (points.size() <= 2) {
    return points;
  }

  std::vector<Point> pts = points;
  std::ranges::sort(pts, [](const Point &a, const Point &b) { return (a.x != b.x) ? (a.x < b.x) : (a.y < b.y); });

  auto [first, last] = std::ranges::unique(pts);
  pts.erase(first, last);

  if (pts.size() <= 2) {
    return pts;
  }

  std::vector<Point> lower;
  std::vector<Point> upper;
  lower.reserve(pts.size());
  upper.reserve(pts.size());

  for (const auto &p : pts) {
    while (lower.size() >= 2 && Cross(lower[lower.size() - 2], lower.back(), p) <= 0) {
      lower.pop_back();
    }
    lower.push_back(p);
  }

  for (const auto &p : std::ranges::reverse_view(pts)) {
    while (upper.size() >= 2 && Cross(upper[upper.size() - 2], upper.back(), p) <= 0) {
      upper.pop_back();
    }
    upper.push_back(p);
  }

  lower.pop_back();
  upper.pop_back();
  lower.insert(lower.end(), upper.begin(), upper.end());
  return lower;
}

}  // namespace shkenev_i_constra_hull_for_binary_image
