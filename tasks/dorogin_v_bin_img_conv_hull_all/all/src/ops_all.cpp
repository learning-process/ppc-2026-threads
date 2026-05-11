#include "dorogin_v_bin_img_conv_hull_all/all/include/ops_all.hpp"

#include <mpi.h>
#include <omp.h>

#include <algorithm>
#include <cstddef>
#include <cstdint>
#include <queue>
#include <ranges>
#include <thread>
#include <utility>
#include <vector>

#include "dorogin_v_bin_img_conv_hull_all/common/include/common.hpp"
#include "oneapi/tbb/parallel_for.h"

namespace dorogin_v_bin_img_conv_hull_all {
namespace {

inline bool InBounds(const int x, const int y, const int w, const int h) {
  return (x >= 0) && (y >= 0) && (x < w) && (y < h);
}

inline std::size_t Idx(const int x, const int y, const int w) {
  return (static_cast<std::size_t>(y) * static_cast<std::size_t>(w)) + static_cast<std::size_t>(x);
}

inline std::int64_t Cross(const Point &o, const Point &a, const Point &b) {
  return (static_cast<std::int64_t>(a.x - o.x) * static_cast<std::int64_t>(b.y - o.y)) -
         (static_cast<std::int64_t>(a.y - o.y) * static_cast<std::int64_t>(b.x - o.x));
}

inline std::vector<Point> ConvexHullMonotonicChain(std::vector<Point> pts) {
  if (pts.empty()) {
    return {};
  }
  std::ranges::sort(pts, [](const Point &a, const Point &b) { return (a.x < b.x) || ((a.x == b.x) && (a.y < b.y)); });
  const auto uniq =
      std::ranges::unique(pts, [](const Point &a, const Point &b) { return (a.x == b.x) && (a.y == b.y); });
  pts.erase(uniq.begin(), pts.end());
  if (pts.size() <= 1) {
    return pts;
  }

  std::vector<Point> lower;
  std::vector<Point> upper;
  lower.reserve(pts.size());
  upper.reserve(pts.size());
  for (const auto &p : pts) {
    while ((lower.size() >= 2) && (Cross(lower[lower.size() - 2], lower[lower.size() - 1], p) <= 0)) {
      lower.pop_back();
    }
    lower.push_back(p);
  }
  for (std::size_t i = pts.size(); i-- > 0;) {
    const auto &p = pts[i];
    while ((upper.size() >= 2) && (Cross(upper[upper.size() - 2], upper[upper.size() - 1], p) <= 0)) {
      upper.pop_back();
    }
    upper.push_back(p);
  }

  lower.pop_back();
  upper.pop_back();
  lower.insert(lower.end(), upper.begin(), upper.end());
  return lower;
}

inline void TryPush4(const BinaryImage &img, const int w, const int h, const int nx, const int ny,
                     std::vector<std::uint8_t> &vis, std::queue<Point> &q) {
  if (!InBounds(nx, ny, w, h)) {
    return;
  }
  const std::size_t nid = Idx(nx, ny, w);
  if ((img.data[nid] == 1U) && (vis[nid] == 0U)) {
    vis[nid] = 1U;
    q.push(Point{.x = nx, .y = ny});
  }
}

inline std::vector<Point> BfsComponent4(const BinaryImage &img, const int w, const int h, const int sx, const int sy,
                                        std::vector<std::uint8_t> &vis) {
  std::vector<Point> pts;
  std::queue<Point> q;
  const std::size_t sid = Idx(sx, sy, w);
  vis[sid] = 1U;
  q.push(Point{.x = sx, .y = sy});
  while (!q.empty()) {
    const Point p = q.front();
    q.pop();
    pts.push_back(p);
    TryPush4(img, w, h, p.x + 1, p.y, vis, q);
    TryPush4(img, w, h, p.x - 1, p.y, vis, q);
    TryPush4(img, w, h, p.x, p.y + 1, vis, q);
    TryPush4(img, w, h, p.x, p.y - 1, vis, q);
  }
  return pts;
}

inline std::vector<std::vector<Point>> ExtractComponents4(const BinaryImage &img) {
  std::vector<std::vector<Point>> comps;
  if ((img.width <= 0) || (img.height <= 0)) {
    return comps;
  }
  const std::size_t n = static_cast<std::size_t>(img.width) * static_cast<std::size_t>(img.height);
  std::vector<std::uint8_t> vis(n, 0U);
  for (int row = 0; row < img.height; ++row) {
    for (int col = 0; col < img.width; ++col) {
      const std::size_t id = Idx(col, row, img.width);
      if ((img.data[id] == 0U) || (vis[id] != 0U)) {
        continue;
      }
      comps.push_back(BfsComponent4(img, img.width, img.height, col, row, vis));
    }
  }
  return comps;
}

inline bool ValidateInput(const BinaryImage &img) {
  if ((img.width <= 0) || (img.height <= 0)) {
    return false;
  }
  const std::size_t need = static_cast<std::size_t>(img.width) * static_cast<std::size_t>(img.height);
  return img.data.size() == need;
}

inline BinaryImage NormalizeImage(BinaryImage img) {
  tbb::parallel_for(static_cast<std::size_t>(0), img.data.size(),
                    [&](const std::size_t i) { img.data[i] = (img.data[i] > 0U) ? 1U : 0U; });
  return img;
}

inline void BuildHullsRange(const std::vector<std::vector<Point>> &comps, OutType &hulls, const std::size_t begin,
                            const std::size_t end) {
  for (std::size_t i = begin; i < end; ++i) {
    hulls[i] = ConvexHullMonotonicChain(comps[i]);
  }
}

inline OutType SolveOnRoot(const BinaryImage &img) {
  const BinaryImage normalized = NormalizeImage(img);
  const auto comps = ExtractComponents4(normalized);
  OutType hulls(comps.size());

  if (comps.size() >= 2) {
    const std::size_t mid = comps.size() / 2;
    std::thread left([&]() { BuildHullsRange(comps, hulls, 0, mid); });
    std::thread right([&]() { BuildHullsRange(comps, hulls, mid, comps.size()); });
    left.join();
    right.join();
  } else if (comps.size() == 1) {
    hulls[0] = ConvexHullMonotonicChain(comps[0]);
  }

  const int count = static_cast<int>(hulls.size());
#pragma omp parallel for default(none) shared(hulls) schedule(static)
  for (int i = 0; i < count; ++i) {
    std::ranges::sort(hulls[static_cast<std::size_t>(i)],
                      [](const Point &a, const Point &b) { return (a.x < b.x) || ((a.x == b.x) && (a.y < b.y)); });
  }

  return hulls;
}

inline void PackHulls(const OutType &hulls, std::vector<int> &counts, std::vector<int> &coords) {
  counts.clear();
  coords.clear();
  counts.reserve(hulls.size());
  for (const auto &hull : hulls) {
    counts.push_back(static_cast<int>(hull.size()));
    for (const auto &point : hull) {
      coords.push_back(point.x);
      coords.push_back(point.y);
    }
  }
}

inline OutType UnpackHulls(const int hull_count, const std::vector<int> &counts, const std::vector<int> &coords) {
  OutType hulls;
  hulls.reserve(static_cast<std::size_t>(std::max(0, hull_count)));
  std::size_t coord_pos = 0;
  for (int hull_idx = 0; hull_idx < hull_count; ++hull_idx) {
    const int point_count = counts[static_cast<std::size_t>(hull_idx)];
    std::vector<Point> hull;
    hull.reserve(static_cast<std::size_t>(std::max(0, point_count)));
    for (int point_idx = 0; point_idx < point_count; ++point_idx) {
      hull.emplace_back(coords[coord_pos], coords[coord_pos + 1]);
      coord_pos += 2;
    }
    hulls.push_back(std::move(hull));
  }
  return hulls;
}

inline void BroadcastHulls(const int rank, OutType &hulls) {
  int hull_count = 0;
  if (rank == 0) {
    hull_count = static_cast<int>(hulls.size());
  }
  MPI_Bcast(&hull_count, 1, MPI_INT, 0, MPI_COMM_WORLD);

  std::vector<int> counts;
  std::vector<int> coords;
  if (rank == 0) {
    PackHulls(hulls, counts, coords);
  } else {
    counts.resize(static_cast<std::size_t>(std::max(0, hull_count)));
  }

  if (hull_count > 0) {
    MPI_Bcast(counts.data(), hull_count, MPI_INT, 0, MPI_COMM_WORLD);
  }

  int coord_count = 0;
  if (rank == 0) {
    coord_count = static_cast<int>(coords.size());
  }
  MPI_Bcast(&coord_count, 1, MPI_INT, 0, MPI_COMM_WORLD);
  if (rank != 0) {
    coords.resize(static_cast<std::size_t>(coord_count));
  }
  if (coord_count > 0) {
    MPI_Bcast(coords.data(), coord_count, MPI_INT, 0, MPI_COMM_WORLD);
  }

  hulls = UnpackHulls(hull_count, counts, coords);
}

}  // namespace

DoroginVConvHullAll::DoroginVConvHullAll(const InType &in) {
  SetTypeOfTask(GetStaticTypeOfTask());
  GetInput() = in;
  GetOutput().clear();
}

bool DoroginVConvHullAll::ValidationImpl() {
  return ValidateInput(GetInput());
}

bool DoroginVConvHullAll::PreProcessingImpl() {
  local_out_.clear();
  return true;
}

bool DoroginVConvHullAll::RunImpl() {
  int rank = 0;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  int is_valid = 0;
  if (rank == 0) {
    is_valid = ValidationImpl() ? 1 : 0;
  }
  MPI_Bcast(&is_valid, 1, MPI_INT, 0, MPI_COMM_WORLD);
  if (is_valid == 0) {
    return false;
  }

  if (rank == 0) {
    local_out_ = SolveOnRoot(GetInput());
  } else {
    local_out_.clear();
  }

  BroadcastHulls(rank, local_out_);
  MPI_Barrier(MPI_COMM_WORLD);
  return true;
}

bool DoroginVConvHullAll::PostProcessingImpl() {
  GetOutput() = local_out_;
  return true;
}

}  // namespace dorogin_v_bin_img_conv_hull_all
