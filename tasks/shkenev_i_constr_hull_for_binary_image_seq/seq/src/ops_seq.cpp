#include "shkenev_i_constr_hull_for_binary_image_seq/seq/include/ops_seq.hpp"

#include <algorithm>
#include <array>
#include <cstddef>
#include <cstdint>
#include <queue>
#include <ranges>
#include <utility>
#include <vector>

#include "shkenev_i_constr_hull_for_binary_image_seq/common/include/common.hpp"

namespace shkenev_i_constr_hull_for_binary_image_seq {

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

ShkenevIConstrHullSeq::ShkenevIConstrHullSeq(const InType &in) : work_(in) {
  SetTypeOfTask(GetStaticTypeOfTask());
  GetInput() = in;
}

bool ShkenevIConstrHullSeq::ValidationImpl() {
  const auto &in = GetInput();
  bool dims = in.width > 0 && in.height > 0;
  bool size = in.pixels.size() == static_cast<size_t>(in.width) * static_cast<size_t>(in.height);
  return dims && size;
}

bool ShkenevIConstrHullSeq::PreProcessingImpl() {
  work_ = GetInput();
  ThresholdImage();
  return true;
}

bool ShkenevIConstrHullSeq::RunImpl() {
  FindComponents();

  work_.convex_hulls.clear();
  work_.convex_hulls.reserve(work_.components.size());

  for (const auto &comp : work_.components) {
    if (comp.empty()) {
      continue;
    }
    if (comp.size() <= 2) {
      work_.convex_hulls.push_back(comp);
    } else {
      work_.convex_hulls.push_back(BuildHull(comp));
    }
  }

  GetOutput() = work_;
  return true;
}

bool ShkenevIConstrHullSeq::PostProcessingImpl() {
  return true;
}

size_t ShkenevIConstrHullSeq::Index(int x, int y, int width) {
  return (static_cast<size_t>(y) * static_cast<size_t>(width)) + static_cast<size_t>(x);
}

void ShkenevIConstrHullSeq::ThresholdImage() {
  for (auto &p : work_.pixels) {
    p = IsForeground(p) ? static_cast<uint8_t>(255) : static_cast<uint8_t>(0);
  }
}

void ShkenevIConstrHullSeq::ExploreComponent(int start_col, int start_row, int width, int height,
                                             std::vector<bool> &visited, std::vector<Point> &component) {
  std::queue<Point> queue;
  queue.emplace(start_col, start_row);
  visited[Index(start_col, start_row, width)] = true;

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
      if (visited[next_idx] || work_.pixels[next_idx] == 0) {
        continue;
      }

      visited[next_idx] = true;
      queue.emplace(next_x, next_y);
    }
  }
}

void ShkenevIConstrHullSeq::FindComponents() {
  int width = work_.width;
  int height = work_.height;

  std::vector<bool> visited(static_cast<size_t>(width) * static_cast<size_t>(height), false);
  work_.components.clear();

  for (int row = 0; row < height; ++row) {
    for (int col = 0; col < width; ++col) {
      size_t idx = Index(col, row, width);
      if (work_.pixels[idx] == 0 || visited[idx]) {
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

std::vector<Point> ShkenevIConstrHullSeq::BuildHull(const std::vector<Point> &points) {
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

}  // namespace shkenev_i_constr_hull_for_binary_image_seq
