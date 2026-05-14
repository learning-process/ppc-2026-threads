#include "dorogin_v_bin_img_conv_hull/stl/include/ops_stl.hpp"

#include <algorithm>
#include <array>
#include <cstddef>
#include <cstdint>
#include <functional>
#include <ranges>
#include <stack>
#include <thread>
#include <utility>
#include <vector>

#include "dorogin_v_bin_img_conv_hull/common/include/common.hpp"

namespace dorogin_v_bin_img_conv_hull {

namespace {

constexpr std::uint8_t kThreshold = 128;

inline bool IsInside(int col, int row, int width, int height) {
  return col >= 0 && row >= 0 && col < width && row < height;
}

std::int64_t Cross(const Point &a, const Point &b, const Point &c) {
  const std::int64_t x1 = static_cast<std::int64_t>(b.x) - static_cast<std::int64_t>(a.x);
  const std::int64_t y1 = static_cast<std::int64_t>(b.y) - static_cast<std::int64_t>(a.y);
  const std::int64_t x2 = static_cast<std::int64_t>(c.x) - static_cast<std::int64_t>(a.x);
  const std::int64_t y2 = static_cast<std::int64_t>(c.y) - static_cast<std::int64_t>(a.y);
  return (x1 * y2) - (y1 * x2);
}

}  // namespace

DoroginVBinImgConvHullSTL::DoroginVBinImgConvHullSTL(const InType &in) : w_(in) {
  SetTypeOfTask(GetStaticTypeOfTask());
  GetInput() = in;
}

bool DoroginVBinImgConvHullSTL::ValidationImpl() {
  const auto &in = GetInput();
  if (in.width <= 0 || in.height <= 0) {
    return false;
  }
  return in.pixels.size() == static_cast<std::size_t>(in.width) * static_cast<std::size_t>(in.height);
}

bool DoroginVBinImgConvHullSTL::PreProcessingImpl() {
  w_ = GetInput();
  ThresholdImage();
  return true;
}

void DoroginVBinImgConvHullSTL::FillConvexHullsChunk(const std::vector<std::vector<Point>> &components,
                                                     std::vector<std::vector<Point>> &convex_hulls,
                                                     const std::size_t begin, const std::size_t end) {
  for (std::size_t i = begin; i < end; ++i) {
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
}

bool DoroginVBinImgConvHullSTL::RunImpl() {
  FindComponents();

  w_.convex_hulls.clear();
  w_.convex_hulls.resize(w_.components.size());

  auto &components = w_.components;
  auto &convex_hulls = w_.convex_hulls;

  const std::size_t n = components.size();
  const std::size_t num_threads =
      std::max(std::size_t{1}, static_cast<std::size_t>(std::thread::hardware_concurrency()));
  const std::size_t chunk = (n + num_threads - 1) / num_threads;

  std::vector<std::thread> threads;
  threads.reserve(num_threads);

  for (std::size_t thread_index = 0; thread_index < num_threads; ++thread_index) {
    const std::size_t begin = thread_index * chunk;
    if (begin >= n) {
      break;
    }
    const std::size_t end = std::min(begin + chunk, n);
    threads.emplace_back([&, begin, end]() { FillConvexHullsChunk(components, convex_hulls, begin, end); });
  }

  for (auto &worker : threads) {
    worker.join();
  }

  GetOutput() = w_;
  return true;
}

bool DoroginVBinImgConvHullSTL::PostProcessingImpl() {
  return true;
}

std::size_t DoroginVBinImgConvHullSTL::Index(int col, int row, int width) {
  return (static_cast<std::size_t>(row) * static_cast<std::size_t>(width)) + static_cast<std::size_t>(col);
}

void DoroginVBinImgConvHullSTL::ThresholdImage() {
  std::ranges::transform(w_.pixels, w_.pixels.begin(), [](std::uint8_t p) {
    return p > kThreshold ? static_cast<std::uint8_t>(255) : static_cast<std::uint8_t>(0);
  });
}

void DoroginVBinImgConvHullSTL::ExploreComponent(int start_col, int start_row, int width, int height,
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

    for (std::size_t dir = 0; dir < dx.size(); ++dir) {
      const int nx = p.x + dx.at(dir);
      const int ny = p.y + dy.at(dir);

      if (!IsInside(nx, ny, width, height)) {
        continue;
      }

      const std::size_t idx = Index(nx, ny, width);
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

void DoroginVBinImgConvHullSTL::FindComponents() {
  const int width = w_.width;
  const int height = w_.height;

  std::vector<bool> visited(static_cast<std::size_t>(width) * static_cast<std::size_t>(height), false);
  w_.components.clear();

  for (int row = 0; row < height; ++row) {
    for (int col = 0; col < width; ++col) {
      const std::size_t idx = Index(col, row, width);

      if (w_.pixels[idx] == 0 || visited[idx]) {
        continue;
      }

      std::vector<Point> comp;
      ExploreComponent(col, row, width, height, visited, comp);

      w_.components.push_back(std::move(comp));
    }
  }
}

std::vector<Point> DoroginVBinImgConvHullSTL::BuildHull(const std::vector<Point> &points) {
  std::vector<Point> pts = points;

  std::ranges::sort(pts, std::less<>{});
  const auto uniq = std::ranges::unique(pts);
  pts.erase(uniq.begin(), uniq.end());

  if (pts.size() <= 2) {
    return pts;
  }

  std::vector<Point> hull;

  for (const auto &p : pts) {
    while (hull.size() >= 2 && Cross(hull[hull.size() - 2], hull.back(), p) <= 0) {
      hull.pop_back();
    }
    hull.push_back(p);
  }

  const std::size_t lsize = hull.size();
  for (int i = static_cast<int>(pts.size()) - 2; i >= 0; --i) {
    while (hull.size() > lsize && Cross(hull[hull.size() - 2], hull.back(), pts[static_cast<std::size_t>(i)]) <= 0) {
      hull.pop_back();
    }
    hull.push_back(pts[static_cast<std::size_t>(i)]);
  }

  hull.pop_back();
  return hull;
}

}  // namespace dorogin_v_bin_img_conv_hull
