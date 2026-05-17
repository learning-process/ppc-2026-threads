#include "shkenev_i_constra_hull_for_binary_image/all/include/ops_all.hpp"

#include <mpi.h>
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

constexpr std::uint8_t kThreshold = 128;
constexpr std::array<std::pair<int, int>, 4> kDirs = {{{1, 0}, {-1, 0}, {0, 1}, {0, -1}}};

bool InBounds(int x, int y, int width, int height) {
  return x >= 0 && x < width && y >= 0 && y < height;
}

bool IsForeground(std::uint8_t pixel) {
  return pixel > kThreshold;
}

std::int64_t Cross(const Point &a, const Point &b, const Point &c) {
  return (static_cast<std::int64_t>(b.x - a.x) * static_cast<std::int64_t>(c.y - b.y)) -
         (static_cast<std::int64_t>(b.y - a.y) * static_cast<std::int64_t>(c.x - b.x));
}

std::vector<int> SerializeHulls(const std::vector<std::vector<Point>> &hulls) {
  std::vector<int> flat;
  flat.reserve(1);
  flat.push_back(static_cast<int>(hulls.size()));
  for (const auto &hull : hulls) {
    flat.push_back(static_cast<int>(hull.size()));
    for (const auto &point : hull) {
      flat.push_back(point.x);
      flat.push_back(point.y);
    }
  }
  return flat;
}

std::vector<std::vector<Point>> DeserializeHulls(const std::vector<int> &flat) {
  if (flat.empty()) {
    return {};
  }

  std::size_t pos = 0;
  const int hull_count = flat[pos++];
  std::vector<std::vector<Point>> hulls;
  hulls.reserve(static_cast<std::size_t>(std::max(0, hull_count)));

  for (int hull_idx = 0; hull_idx < hull_count && pos < flat.size(); ++hull_idx) {
    const int point_count = flat[pos++];
    std::vector<Point> hull;
    hull.reserve(static_cast<std::size_t>(std::max(0, point_count)));

    for (int point_idx = 0; point_idx < point_count && (pos + 1) < flat.size(); ++point_idx) {
      hull.emplace_back(flat[pos], flat[pos + 1]);
      pos += 2;
    }
    hulls.push_back(std::move(hull));
  }

  return hulls;
}

}  // namespace

size_t ShkenevIConstrHullALL::Index(int x, int y, int width) {
  return (static_cast<size_t>(y) * static_cast<size_t>(width)) + static_cast<size_t>(x);
}

ShkenevIConstrHullALL::ShkenevIConstrHullALL(const InType &in) : work_(in) {
  SetTypeOfTask(GetStaticTypeOfTask());
  GetInput() = in;
}

bool ShkenevIConstrHullALL::ValidationImpl() {
  const auto &in = GetInput();
  return in.width > 0 && in.height > 0 &&
         in.pixels.size() == static_cast<size_t>(in.width) * static_cast<size_t>(in.height);
}

bool ShkenevIConstrHullALL::PreProcessingImpl() {
  work_ = GetInput();
  work_.components.clear();
  work_.convex_hulls.clear();

  auto &pixels = work_.pixels;
  tbb::parallel_for(static_cast<size_t>(0), pixels.size(), [&](size_t idx) {
    pixels[idx] = IsForeground(pixels[idx]) ? static_cast<std::uint8_t>(255) : static_cast<std::uint8_t>(0);
  });

  return true;
}

void ShkenevIConstrHullALL::BFS(int sx, int sy, int width, int height, std::vector<std::uint8_t> &visited,
                                std::vector<Point> &component) {
  std::queue<Point> queue;
  queue.emplace(sx, sy);
  visited[Index(sx, sy, width)] = 1;

  while (!queue.empty()) {
    const Point current = queue.front();
    queue.pop();
    component.push_back(current);

    for (auto [dx, dy] : kDirs) {
      const int next_x = current.x + dx;
      const int next_y = current.y + dy;
      if (!InBounds(next_x, next_y, width, height)) {
        continue;
      }

      const size_t next_idx = Index(next_x, next_y, width);
      if (visited[next_idx] != 0U || work_.pixels[next_idx] == 0) {
        continue;
      }

      visited[next_idx] = 1;
      queue.emplace(next_x, next_y);
    }
  }
}

void ShkenevIConstrHullALL::FindComponentsTBB() {
  const int width = work_.width;
  const int height = work_.height;

  std::vector<std::uint8_t> visited(static_cast<size_t>(width) * static_cast<size_t>(height), 0);
  work_.components.clear();

  for (int row = 0; row < height; ++row) {
    for (int col = 0; col < width; ++col) {
      const size_t idx = Index(col, row, width);
      if (work_.pixels[idx] == 0 || visited[idx] != 0U) {
        continue;
      }

      std::vector<Point> component;
      BFS(col, row, width, height, visited, component);
      if (!component.empty()) {
        work_.components.push_back(std::move(component));
      }
    }
  }
}

std::vector<Point> ShkenevIConstrHullALL::BuildHull(const std::vector<Point> &points) {
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

  for (const auto &point : pts) {
    while (lower.size() >= 2 && Cross(lower[lower.size() - 2], lower.back(), point) <= 0) {
      lower.pop_back();
    }
    lower.push_back(point);
  }

  for (const auto &point : std::ranges::reverse_view(pts)) {
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

bool ShkenevIConstrHullALL::RunImpl() {
  int rank = 0;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  if (rank == 0) {
    FindComponentsTBB();

    auto &components = work_.components;
    auto &hulls = work_.convex_hulls;
    hulls.resize(components.size());

    tbb::parallel_for(static_cast<size_t>(0), components.size(), [&](size_t idx) {
      const auto &component = components[idx];
      hulls[idx] = (component.size() <= 2) ? component : BuildHull(component);
    });
  }

  std::vector<int> flat_hulls;
  if (rank == 0) {
    flat_hulls = SerializeHulls(work_.convex_hulls);
  }

  int flat_size = static_cast<int>(flat_hulls.size());
  MPI_Bcast(&flat_size, 1, MPI_INT, 0, MPI_COMM_WORLD);
  if (rank != 0) {
    flat_hulls.resize(static_cast<std::size_t>(flat_size));
  }
  MPI_Bcast(flat_hulls.data(), flat_size, MPI_INT, 0, MPI_COMM_WORLD);

  work_.convex_hulls = DeserializeHulls(flat_hulls);
  work_.components.clear();
  MPI_Barrier(MPI_COMM_WORLD);
  GetOutput() = work_;
  return true;
}

bool ShkenevIConstrHullALL::PostProcessingImpl() {
  return true;
}

}  // namespace shkenev_i_constra_hull_for_binary_image
