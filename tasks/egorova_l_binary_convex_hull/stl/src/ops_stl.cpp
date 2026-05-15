#include "egorova_l_binary_convex_hull/stl/include/ops_stl.hpp"

#include <algorithm>
#include <atomic>
#include <cstddef>
#include <cstdint>
#include <queue>
#include <thread>
#include <tuple>
#include <utility>
#include <vector>

#include "egorova_l_binary_convex_hull/common/include/common.hpp"
#include "util/include/util.hpp"

namespace egorova_l_binary_convex_hull {

namespace {

std::vector<Point> RunBFS(int start_col, int start_row, int w, int h, const std::vector<uint8_t> &data,
                          std::vector<bool> &visited) {
  static const std::vector<std::pair<int, int>> kDirs = {{0, 1}, {0, -1}, {1, 0}, {-1, 0}};

  std::vector<Point> component;
  std::queue<std::pair<int, int>> q;

  q.emplace(start_col, start_row);
  visited[(static_cast<size_t>(start_row) * w) + start_col] = true;

  while (!q.empty()) {
    auto [cx, cy] = q.front();
    q.pop();
    component.push_back({cx, cy});

    for (const auto &dir : kDirs) {
      int nx = cx + dir.first;
      int ny = cy + dir.second;

      if (nx < 0 || nx >= w || ny < 0 || ny >= h) {
        continue;
      }

      const size_t nidx = (static_cast<size_t>(ny) * w) + nx;
      if (data[nidx] != 0 && !visited[nidx]) {
        visited[nidx] = true;
        q.emplace(nx, ny);
      }
    }
  }
  return component;
}

int64_t Cross(const Point &o, const Point &a, const Point &b) {
  return (static_cast<int64_t>(a.x - o.x) * static_cast<int64_t>(b.y - o.y)) -
         (static_cast<int64_t>(a.y - o.y) * static_cast<int64_t>(b.x - o.x));
}

}  // namespace

BinaryConvexHullSTL::BinaryConvexHullSTL(const InType &in) {
  SetTypeOfTask(GetStaticTypeOfTask());
  GetInput() = in;
  GetOutput() = {};
}

bool BinaryConvexHullSTL::ValidationImpl() {
  const auto &img = GetInput();
  if (img.width <= 0 || img.height <= 0) {
    return false;
  }
  return static_cast<int>(img.data.size()) == img.width * img.height;
}

bool BinaryConvexHullSTL::PreProcessingImpl() {
  components_ = FindComponents(GetInput());
  return true;
}

bool BinaryConvexHullSTL::RunImpl() {
  const int n = static_cast<int>(components_.size());
  if (n == 0) {
    GetOutput().clear();
    return true;
  }

  const int num_threads = std::min(ppc::util::GetNumThreads(), n);
  std::vector<std::vector<Point>> result(n);
  std::atomic<int> next_idx(0);

  auto worker = [&]() {
    int i = 0;
    while ((i = next_idx.fetch_add(1)) < n) {
      result[i] = ConvexHull(components_[i]);
    }
  };

  std::vector<std::thread> threads;
  threads.reserve(num_threads);
  for (int i = 0; i < num_threads; ++i) {
    threads.emplace_back(worker);
  }

  for (auto &t : threads) {
    if (t.joinable()) {
      t.join();
    }
  }

  GetOutput() = std::move(result);
  return true;
}

bool BinaryConvexHullSTL::PostProcessingImpl() {
  return true;
}

std::vector<std::vector<Point>> BinaryConvexHullSTL::FindComponents(const ImageData &img) {
  const int w = img.width;
  const int h = img.height;
  std::vector<bool> visited(static_cast<size_t>(w) * h, false);
  std::vector<std::vector<Point>> components;

  for (int row = 0; row < h; ++row) {
    for (int col = 0; col < w; ++col) {
      const size_t idx = (static_cast<size_t>(row) * w) + col;
      if (img.data[idx] == 0 || visited[idx]) {
        continue;
      }
      components.push_back(RunBFS(col, row, w, h, img.data, visited));
    }
  }

  return components;
}

std::vector<Point> BinaryConvexHullSTL::ConvexHull(std::vector<Point> points) {
  const int n = static_cast<int>(points.size());
  if (n == 0) {
    return {};
  }
  if (n == 1) {
    return points;
  }

  std::ranges::sort(points, [](const Point &a, const Point &b) { return std::tie(a.x, a.y) < std::tie(b.x, b.y); });

  auto [first, last] =
      std::ranges::unique(points, [](const Point &a, const Point &b) { return a.x == b.x && a.y == b.y; });
  points.erase(first, last);

  const int m = static_cast<int>(points.size());
  if (m <= 2) {
    return points;
  }

  std::vector<Point> hull;
  hull.reserve(static_cast<size_t>(2) * static_cast<size_t>(m));

  for (int i = 0; i < m; ++i) {
    while (hull.size() >= 2 && Cross(hull[hull.size() - 2], hull.back(), points[i]) <= 0) {
      hull.pop_back();
    }
    hull.push_back(points[i]);
  }

  const size_t lower_hull_size = hull.size();
  for (int i = m - 2; i >= 0; --i) {
    while (hull.size() > lower_hull_size && Cross(hull[hull.size() - 2], hull.back(), points[i]) <= 0) {
      hull.pop_back();
    }
    hull.push_back(points[i]);
  }
  hull.pop_back();

  return hull;
}

}  // namespace egorova_l_binary_convex_hull
