#include "peryashkin_v_binary_component_contour_processing/all/include/ops_all.hpp"

#include <mpi.h>
#include <omp.h>

#include <algorithm>
#include <cstddef>
#include <cstdint>
#include <queue>
#include <ranges>
#include <utility>
#include <vector>

#include "peryashkin_v_binary_component_contour_processing/common/include/common.hpp"

namespace peryashkin_v_binary_component_contour_processing {

namespace {

inline bool InBounds(int x, int y, int w, int h) {
  return (x >= 0) && (y >= 0) && (x < w) && (y < h);
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

  if (pts.size() == 1) {
    return pts;
  }

  std::vector<Point> lower;
  lower.reserve(pts.size());
  for (const auto &p : pts) {
    while ((lower.size() >= 2) && (Cross(lower[lower.size() - 2], lower[lower.size() - 1], p) <= 0)) {
      lower.pop_back();
    }
    lower.push_back(p);
  }

  std::vector<Point> upper;
  upper.reserve(pts.size());
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

inline std::size_t Idx(int x, int y, int w) {
  return (static_cast<std::size_t>(y) * static_cast<std::size_t>(w)) + static_cast<std::size_t>(x);
}

inline void TryPush4(const BinaryImage &img, int w, int h, int nx, int ny, std::vector<std::uint8_t> &vis,
                     std::queue<Point> &q) {
  if (!InBounds(nx, ny, w, h)) {
    return;
  }
  const std::size_t nid = Idx(nx, ny, w);
  if ((img.data[nid] == 1) && (vis[nid] == 0U)) {
    vis[nid] = 1U;
    q.push(Point{.x = nx, .y = ny});
  }
}

inline std::vector<Point> BfsComponent4(const BinaryImage &img, int w, int h, int sx, int sy,
                                        std::vector<std::uint8_t> &vis) {
  std::vector<Point> pts;
  pts.reserve(128);

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
  const int w = img.width;
  const int h = img.height;

  std::vector<std::vector<Point>> comps;
  if ((w <= 0) || (h <= 0)) {
    return comps;
  }

  const std::size_t n = static_cast<std::size_t>(w) * static_cast<std::size_t>(h);
  std::vector<std::uint8_t> vis(n, 0U);

  for (int yy = 0; yy < h; ++yy) {
    for (int xx = 0; xx < w; ++xx) {
      const std::size_t id = Idx(xx, yy, w);
      if ((img.data[id] == 0) || (vis[id] != 0U)) {
        continue;
      }
      comps.push_back(BfsComponent4(img, w, h, xx, yy, vis));
    }
  }

  return comps;
}

inline std::vector<int> FlattenComponents(const std::vector<std::vector<Point>> &comps) {
  std::vector<int> flat;
  for (const auto &comp : comps) {
    flat.push_back(static_cast<int>(comp.size()));
    for (const auto &p : comp) {
      flat.push_back(p.x);
      flat.push_back(p.y);
    }
  }
  return flat;
}

inline std::vector<std::vector<Point>> UnflattenComponents(const std::vector<int> &flat) {
  std::vector<std::vector<Point>> comps;
  std::size_t pos = 0;

  while (pos < flat.size()) {
    const int cnt = flat[pos++];
    std::vector<Point> comp;
    comp.reserve(static_cast<std::size_t>(cnt));

    for (int i = 0; i < cnt; ++i) {
      const int x = flat[pos++];
      const int y = flat[pos++];
      comp.push_back(Point{.x = x, .y = y});
    }

    comps.push_back(std::move(comp));
  }

  return comps;
}

inline std::vector<int> MakeDispls(const std::vector<int> &counts) {
  std::vector<int> displs(counts.size(), 0);
  for (std::size_t i = 1; i < counts.size(); ++i) {
    displs[i] = displs[i - 1] + counts[i - 1];
  }
  return displs;
}

inline std::vector<int> MakeComponentCounts(int total_components, int proc_count) {
  std::vector<int> comp_counts(proc_count, 0);
  for (int i = 0; i < proc_count; ++i) {
    comp_counts[i] = (total_components / proc_count) + ((i < (total_components % proc_count)) ? 1 : 0);
  }
  return comp_counts;
}

inline std::vector<int> FlattenDistributedComponents(const std::vector<std::vector<Point>> &all_components,
                                                     const std::vector<int> &comp_counts,
                                                     std::vector<int> &send_counts) {
  std::vector<int> flat_send;
  int comp_offset = 0;

  for (std::size_t proc = 0; proc < comp_counts.size(); ++proc) {
    std::vector<std::vector<Point>> part;
    part.reserve(static_cast<std::size_t>(comp_counts[proc]));

    for (int j = 0; j < comp_counts[proc]; ++j) {
      const auto comp_index = static_cast<std::size_t>(comp_offset) + static_cast<std::size_t>(j);
      part.push_back(all_components[comp_index]);
    }

    comp_offset += comp_counts[proc];

    std::vector<int> flat_part = FlattenComponents(part);
    send_counts[proc] = static_cast<int>(flat_part.size());
    flat_send.insert(flat_send.end(), flat_part.begin(), flat_part.end());
  }

  return flat_send;
}

inline OutType BuildLocalHulls(std::vector<std::vector<Point>> local_components) {
  OutType local_hulls(local_components.size());
  const int local_components_count = static_cast<int>(local_components.size());

#pragma omp parallel for default(none) shared(local_components, local_hulls, local_components_count)
  for (int i = 0; i < local_components_count; ++i) {
    const auto idx = static_cast<std::size_t>(i);
    local_hulls[idx] = ConvexHullMonotonicChain(std::move(local_components[idx]));
  }

  return local_hulls;
}

inline std::vector<int> GatherFlatHulls(const std::vector<int> &flat_local_hulls, int rank, int size) {
  const int local_hulls_size = static_cast<int>(flat_local_hulls.size());

  std::vector<int> recv_hull_counts(static_cast<std::size_t>(size), 0);
  MPI_Gather(&local_hulls_size, 1, MPI_INT, recv_hull_counts.data(), 1, MPI_INT, 0, MPI_COMM_WORLD);

  std::vector<int> recv_hull_displs;
  std::vector<int> flat_hulls_root;

  if (rank == 0) {
    recv_hull_displs = MakeDispls(recv_hull_counts);

    int total_size = 0;
    if (!recv_hull_counts.empty()) {
      const auto last_idx = recv_hull_counts.size() - 1;
      total_size = recv_hull_displs[last_idx] + recv_hull_counts[last_idx];
    }

    flat_hulls_root.resize(static_cast<std::size_t>(total_size));
  }

  MPI_Gatherv(flat_local_hulls.data(), local_hulls_size, MPI_INT, rank == 0 ? flat_hulls_root.data() : nullptr,
              recv_hull_counts.data(), rank == 0 ? recv_hull_displs.data() : nullptr, MPI_INT, 0, MPI_COMM_WORLD);

  return flat_hulls_root;
}

// NOLINTNEXTLINE(readability-function-cognitive-complexity)
inline OutType SolveALL(const BinaryImage &img) {
  int rank = 0;
  int size = 1;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &size);

  std::vector<int> send_counts(static_cast<std::size_t>(size), 0);
  std::vector<int> send_displs(static_cast<std::size_t>(size), 0);
  std::vector<int> flat_send;

  if (rank == 0) {
    const auto all_components = ExtractComponents4(img);
    const int total_components = static_cast<int>(all_components.size());
    const auto comp_counts = MakeComponentCounts(total_components, size);

    flat_send = FlattenDistributedComponents(all_components, comp_counts, send_counts);
    send_displs = MakeDispls(send_counts);
  }

  int recv_count = 0;
  MPI_Scatter(send_counts.data(), 1, MPI_INT, &recv_count, 1, MPI_INT, 0, MPI_COMM_WORLD);

  std::vector<int> flat_recv(static_cast<std::size_t>(recv_count));
  MPI_Scatterv(rank == 0 ? flat_send.data() : nullptr, send_counts.data(), send_displs.data(), MPI_INT,
               flat_recv.data(), recv_count, MPI_INT, 0, MPI_COMM_WORLD);

  auto local_components = UnflattenComponents(flat_recv);
  const auto local_hulls = BuildLocalHulls(std::move(local_components));
  const auto flat_local_hulls = FlattenComponents(local_hulls);

  std::vector<int> flat_broadcast;
  if (rank == 0) {
    flat_broadcast = GatherFlatHulls(flat_local_hulls, rank, size);
  } else {
    GatherFlatHulls(flat_local_hulls, rank, size);
  }

  int broadcast_size = 0;
  if (rank == 0) {
    broadcast_size = static_cast<int>(flat_broadcast.size());
  }

  MPI_Bcast(&broadcast_size, 1, MPI_INT, 0, MPI_COMM_WORLD);

  if (rank != 0) {
    flat_broadcast.resize(static_cast<std::size_t>(broadcast_size));
  }

  MPI_Bcast(flat_broadcast.data(), broadcast_size, MPI_INT, 0, MPI_COMM_WORLD);

  return UnflattenComponents(flat_broadcast);
}

}  // namespace

PeryashkinVBinaryComponentContourProcessingALL::PeryashkinVBinaryComponentContourProcessingALL(const InType &in) {
  SetTypeOfTask(GetStaticTypeOfTask());
  GetInput() = in;
  GetOutput().clear();
}

bool PeryashkinVBinaryComponentContourProcessingALL::ValidationImpl() {
  const auto &in = GetInput();
  if ((in.width <= 0) || (in.height <= 0)) {
    return false;
  }

  const std::size_t need = static_cast<std::size_t>(in.width) * static_cast<std::size_t>(in.height);
  return in.data.size() == need;
}

bool PeryashkinVBinaryComponentContourProcessingALL::PreProcessingImpl() {
  local_out_.clear();
  return true;
}

bool PeryashkinVBinaryComponentContourProcessingALL::RunImpl() {
  if (!ValidationImpl()) {
    return false;
  }

  local_out_ = SolveALL(GetInput());
  return true;
}

bool PeryashkinVBinaryComponentContourProcessingALL::PostProcessingImpl() {
  GetOutput() = local_out_;
  return true;
}

}  // namespace peryashkin_v_binary_component_contour_processing
