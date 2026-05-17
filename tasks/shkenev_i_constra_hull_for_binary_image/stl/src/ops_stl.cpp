#include "shkenev_i_constra_hull_for_binary_image/stl/include/ops_stl.hpp"

#include <algorithm>
#include <array>
#include <atomic>
#include <cstddef>
#include <cstdint>
#include <thread>
#include <utility>
#include <vector>

#include "shkenev_i_constra_hull_for_binary_image/common/include/common.hpp"
#include "util/include/util.hpp"

namespace shkenev_i_constra_hull_for_binary_image {

namespace {

// 1
constexpr uint8_t kThreshold = 128;

constexpr std::array<std::pair<int, int>, 4> kDirs = {std::pair{1, 0}, std::pair{-1, 0}, std::pair{0, 1},
                                                      std::pair{0, -1}};

inline bool InBounds(int x, int y, int w, int h) {
  return x >= 0 && x < w && y >= 0 && y < h;
}

inline int64_t Cross(const Point &a, const Point &b, const Point &c) {
  return (static_cast<int64_t>(b.x - a.x) * static_cast<int64_t>(c.y - a.y)) -
         (static_cast<int64_t>(b.y - a.y) * static_cast<int64_t>(c.x - a.x));
}

}  // namespace

ShkenevIConstrHullSTL::ShkenevIConstrHullSTL(const InType &in) : work_(in) {
  SetTypeOfTask(GetStaticTypeOfTask());
  GetInput() = in;
}

bool ShkenevIConstrHullSTL::ValidationImpl() {
  const auto &in = GetInput();
  return in.width > 0 && in.height > 0 &&
         in.pixels.size() == static_cast<size_t>(in.width) * static_cast<size_t>(in.height);
}

bool ShkenevIConstrHullSTL::PreProcessingImpl() {
  work_ = GetInput();
  ThresholdImage();
  return true;
}

void ShkenevIConstrHullSTL::ThresholdImage() {
  for (auto &pixel : work_.pixels) {
    pixel = (pixel > kThreshold) ? 255 : 0;
  }
}

size_t ShkenevIConstrHullSTL::Index(int x, int y, int width) {
  return (static_cast<size_t>(y) * static_cast<size_t>(width)) + static_cast<size_t>(x);
}

void ShkenevIConstrHullSTL::ExploreComponent(int start_x, int start_y, int width, int height,
                                             std::vector<uint8_t> &visited, std::vector<Point> &component) {
  std::vector<Point> stack;
  stack.reserve(256);
  component.reserve(256);

  stack.emplace_back(start_x, start_y);
  visited[Index(start_x, start_y, width)] = 1;

  while (!stack.empty()) {
    Point current = stack.back();
    stack.pop_back();

    component.push_back(current);

    for (const auto &[delta_x, delta_y] : kDirs) {
      int neighbor_x = current.x + delta_x;
      int neighbor_y = current.y + delta_y;

      if (!InBounds(neighbor_x, neighbor_y, width, height)) {
        continue;
      }

      size_t idx = Index(neighbor_x, neighbor_y, width);

      if (visited[idx] != 0U || work_.pixels[idx] == 0) {
        continue;
      }

      visited[idx] = 1;
      stack.emplace_back(neighbor_x, neighbor_y);
    }
  }
}

void ShkenevIConstrHullSTL::FindComponents() {
  int width = work_.width;
  int height = work_.height;

  std::vector<uint8_t> visited(static_cast<size_t>(width) * static_cast<size_t>(height), 0);
  work_.components.clear();

  for (int row = 0; row < height; ++row) {
    for (int col = 0; col < width; ++col) {
      size_t idx = Index(col, row, width);

      if (visited[idx] != 0U || work_.pixels[idx] == 0) {
        continue;
      }

      std::vector<Point> component;
      ExploreComponent(col, row, width, height, visited, component);

      if (!component.empty()) {
        work_.components.emplace_back(std::move(component));
      }
    }
  }
}

bool ShkenevIConstrHullSTL::RunImpl() {
  FindComponents();

  auto &components = work_.components;
  auto &hulls = work_.convex_hulls;

  if (components.empty()) {
    GetOutput() = work_;
    return true;
  }

  hulls.resize(components.size());

  int num_threads = std::min<int>(ppc::util::GetNumThreads(), static_cast<int>(components.size()));
  std::vector<std::thread> threads;
  threads.reserve(static_cast<size_t>(num_threads));

  std::atomic<size_t> current_index{0};

  for (int thread_id = 0; thread_id < num_threads; ++thread_id) {
    threads.emplace_back([&]() {
      while (true) {
        size_t idx = current_index.fetch_add(1, std::memory_order_relaxed);
        if (idx >= components.size()) {
          break;
        }

        const auto &component = components[idx];

        if (component.size() <= 2) {
          hulls[idx] = component;
        } else {
          hulls[idx] = BuildHull(component);
        }
      }
    });
  }

  for (auto &thread : threads) {
    thread.join();
  }

  GetOutput() = work_;
  return true;
}

std::vector<Point> ShkenevIConstrHullSTL::BuildHull(const std::vector<Point> &pts_in) {
  std::vector<Point> pts = pts_in;

  std::ranges::sort(pts, [](const Point &a, const Point &b) { return (a.x < b.x) || (a.x == b.x && a.y < b.y); });

  auto [first, last] =
      std::ranges::unique(pts, [](const Point &a, const Point &b) { return a.x == b.x && a.y == b.y; });
  pts.erase(first, last);

  if (pts.size() <= 2) {
    return pts;
  }

  std::vector<Point> lower;
  std::vector<Point> upper;
  lower.reserve(pts.size());
  upper.reserve(pts.size());

  for (const auto &point : pts) {
    while (lower.size() >= 2 && Cross(lower[lower.size() - 2], lower.back(), point) <= 0) {
      lower.pop_back();
    }
    lower.push_back(point);
  }

  for (int i = static_cast<int>(pts.size()) - 1; i >= 0; --i) {
    const auto &point = pts[i];
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

bool ShkenevIConstrHullSTL::PostProcessingImpl() {
  return true;
}

}  // namespace shkenev_i_constra_hull_for_binary_image
