#include "dorogin_v_bin_img_conv_hull/all/include/ops_all.hpp"

#include <mpi.h>

#include <algorithm>
#include <array>
#include <cstddef>
#include <cstdint>
#include <functional>
#include <queue>
#include <ranges>
#include <utility>
#include <vector>

#include "dorogin_v_bin_img_conv_hull/common/include/common.hpp"
#include "oneapi/tbb/parallel_for.h"

namespace dorogin_v_bin_img_conv_hull {

namespace {

constexpr std::uint8_t kThreshold = 128;

constexpr std::array<std::pair<int, int>, 4> kNeighbours = {{{1, 0}, {-1, 0}, {0, 1}, {0, -1}}};

bool CellInImage(int x, int y, int width, int height) {
  return x >= 0 && x < width && y >= 0 && y < height;
}

bool IsForeground(std::uint8_t pixel) {
  return pixel > kThreshold;
}

std::int64_t Orient(const Point &a, const Point &b, const Point &c) {
  const std::int64_t x1 = static_cast<std::int64_t>(b.x) - static_cast<std::int64_t>(a.x);
  const std::int64_t y1 = static_cast<std::int64_t>(b.y) - static_cast<std::int64_t>(a.y);
  const std::int64_t x2 = static_cast<std::int64_t>(c.x) - static_cast<std::int64_t>(a.x);
  const std::int64_t y2 = static_cast<std::int64_t>(c.y) - static_cast<std::int64_t>(a.y);
  return (x1 * y2) - (y1 * x2);
}

std::vector<int> FlattenHulls(const std::vector<std::vector<Point>> &hulls) {
  std::vector<int> buf;
  buf.push_back(static_cast<int>(hulls.size()));
  for (const auto &hull : hulls) {
    buf.push_back(static_cast<int>(hull.size()));
    for (const Point &p : hull) {
      buf.push_back(p.x);
      buf.push_back(p.y);
    }
  }
  return buf;
}

std::vector<std::vector<Point>> RestoreHulls(const std::vector<int> &buf) {
  if (buf.empty()) {
    return {};
  }
  std::size_t pos = 0;
  const int num_hulls = buf[pos++];
  std::vector<std::vector<Point>> hulls;
  hulls.reserve(static_cast<std::size_t>(std::max(0, num_hulls)));
  for (int hull_index = 0; hull_index < num_hulls; ++hull_index) {
    if (pos >= buf.size()) {
      break;
    }
    const int point_count = buf[pos++];
    std::vector<Point> hull;
    hull.reserve(static_cast<std::size_t>(std::max(0, point_count)));
    for (int point_index = 0; point_index < point_count; ++point_index) {
      if (pos + 1 >= buf.size()) {
        break;
      }
      hull.emplace_back(buf[pos], buf[pos + 1]);
      pos += 2;
    }
    hulls.push_back(std::move(hull));
  }
  return hulls;
}

}  // namespace

std::size_t DoroginVBinImgConvHullALL::Index(int col, int row, int width) {
  return (static_cast<std::size_t>(row) * static_cast<std::size_t>(width)) + static_cast<std::size_t>(col);
}

DoroginVBinImgConvHullALL::DoroginVBinImgConvHullALL(const InType &in) : work_(in) {
  SetTypeOfTask(GetStaticTypeOfTask());
  GetInput() = in;
}

bool DoroginVBinImgConvHullALL::ValidationImpl() {
  const auto &in = GetInput();
  if (in.width <= 0 || in.height <= 0) {
    return false;
  }
  return in.pixels.size() == static_cast<std::size_t>(in.width) * static_cast<std::size_t>(in.height);
}

bool DoroginVBinImgConvHullALL::PreProcessingImpl() {
  work_ = GetInput();
  work_.components.clear();
  work_.convex_hulls.clear();
  ThresholdPixelsParallel();
  return true;
}

void DoroginVBinImgConvHullALL::ThresholdPixelsParallel() {
  auto &pixels = work_.pixels;
  tbb::parallel_for(static_cast<std::size_t>(0), pixels.size(), [&](std::size_t idx) {
    pixels[idx] = IsForeground(pixels[idx]) ? static_cast<std::uint8_t>(255) : static_cast<std::uint8_t>(0);
  });
}

void DoroginVBinImgConvHullALL::FloodFill(int seed_x, int seed_y, int width, int height,
                                          std::vector<std::uint8_t> &visited, std::vector<Point> &component) {
  std::queue<Point> q;
  q.emplace(seed_x, seed_y);
  visited[Index(seed_x, seed_y, width)] = 1;
  while (!q.empty()) {
    const Point cur = q.front();
    q.pop();
    component.push_back(cur);
    for (const auto &[dx, dy] : kNeighbours) {
      const int nx = cur.x + dx;
      const int ny = cur.y + dy;
      if (!CellInImage(nx, ny, width, height)) {
        continue;
      }
      const std::size_t ni = Index(nx, ny, width);
      if (visited[ni] != 0U || work_.pixels[ni] == 0) {
        continue;
      }
      visited[ni] = 1;
      q.emplace(nx, ny);
    }
  }
}

void DoroginVBinImgConvHullALL::CollectComponentsSequential() {
  const int width = work_.width;
  const int height = work_.height;
  std::vector<std::uint8_t> visited(static_cast<std::size_t>(width) * static_cast<std::size_t>(height), 0);
  work_.components.clear();

  for (int row = 0; row < height; ++row) {
    for (int col = 0; col < width; ++col) {
      const std::size_t idx = Index(col, row, width);
      if (work_.pixels[idx] == 0 || visited[idx] != 0U) {
        continue;
      }
      std::vector<Point> comp;
      FloodFill(col, row, width, height, visited, comp);
      if (!comp.empty()) {
        work_.components.push_back(std::move(comp));
      }
    }
  }
}

std::vector<Point> DoroginVBinImgConvHullALL::BuildHull(const std::vector<Point> &points) {
  if (points.size() < 3) {
    return points;
  }
  std::vector<Point> pts = points;
  std::ranges::sort(pts, std::less<>{});
  const auto uniq = std::ranges::unique(pts);
  pts.erase(uniq.begin(), uniq.end());
  if (pts.size() < 3) {
    return pts;
  }

  std::vector<Point> lower;
  std::vector<Point> upper;
  lower.reserve(pts.size());
  upper.reserve(pts.size());

  for (const auto &p : pts) {
    while (lower.size() >= 2 && Orient(lower[lower.size() - 2], lower.back(), p) <= 0) {
      lower.pop_back();
    }
    lower.push_back(p);
  }
  for (const auto &p : std::ranges::reverse_view(pts)) {
    while (upper.size() >= 2 && Orient(upper[upper.size() - 2], upper.back(), p) <= 0) {
      upper.pop_back();
    }
    upper.push_back(p);
  }
  lower.pop_back();
  upper.pop_back();
  lower.insert(lower.end(), upper.begin(), upper.end());
  return lower;
}

bool DoroginVBinImgConvHullALL::RunImpl() {
  int rank = 0;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  int ok = 0;
  if (rank == 0) {
    ok = ValidationImpl() ? 1 : 0;
  }
  MPI_Bcast(&ok, 1, MPI_INT, 0, MPI_COMM_WORLD);
  if (ok == 0) {
    return false;
  }

  if (rank == 0) {
    CollectComponentsSequential();
    auto &comps = work_.components;
    auto &hulls = work_.convex_hulls;
    hulls.resize(comps.size());
    tbb::parallel_for(static_cast<std::size_t>(0), comps.size(), [&](std::size_t idx) {
      const auto &comp = comps[idx];
      hulls[idx] = (comp.size() < 3) ? comp : BuildHull(comp);
    });
  }

  std::vector<int> packed;
  if (rank == 0) {
    packed = FlattenHulls(work_.convex_hulls);
  }

  int packed_len = static_cast<int>(packed.size());
  MPI_Bcast(&packed_len, 1, MPI_INT, 0, MPI_COMM_WORLD);
  if (rank != 0) {
    packed.resize(static_cast<std::size_t>(packed_len));
  }
  if (packed_len > 0) {
    MPI_Bcast(packed.data(), packed_len, MPI_INT, 0, MPI_COMM_WORLD);
  }

  work_.convex_hulls = RestoreHulls(packed);
  work_.components.clear();
  MPI_Barrier(MPI_COMM_WORLD);

  GetOutput() = work_;
  return true;
}

bool DoroginVBinImgConvHullALL::PostProcessingImpl() {
  return true;
}

}  // namespace dorogin_v_bin_img_conv_hull
