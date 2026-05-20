#include "kamalagin_a_binary_image_convex_hull/all/include/ops_all.hpp"

#include <mpi.h>
#include <omp.h>

#include <algorithm>
#include <array>
#include <cstddef>
#include <cstdint>
#include <utility>
#include <vector>

#include "kamalagin_a_binary_image_convex_hull/common/include/common.hpp"

namespace kamalagin_a_binary_image_convex_hull {

namespace {

constexpr std::array<std::pair<int, int>, 4> kFourNeighbors = {{{-1, 0}, {1, 0}, {0, -1}, {0, 1}}};

int64_t Cross(const Point &o, const Point &a, const Point &b) {
  return (static_cast<int64_t>(a.x - o.x) * (b.y - o.y)) - (static_cast<int64_t>(a.y - o.y) * (b.x - o.x));
}

void GrahamHull(std::vector<Point> &pts, Hull &out) {
  out.clear();
  const size_t n = pts.size();
  if (n <= 1) {
    if (n == 1) {
      out.push_back(pts[0]);
    }
    return;
  }

  size_t pivot = 0;
  for (size_t i = 1; i < n; ++i) {
    if (pts[i].y < pts[pivot].y || (pts[i].y == pts[pivot].y && pts[i].x < pts[pivot].x)) {
      pivot = i;
    }
  }
  std::swap(pts[0], pts[pivot]);

  const Point origin = pts[0];
  std::sort(pts.begin() + 1, pts.end(), [&origin](const Point &a, const Point &b) {
    const int64_t cr = Cross(origin, a, b);
    if (cr != 0) {
      return cr > 0;
    }
    const int64_t da = (static_cast<int64_t>(a.x - origin.x) * (a.x - origin.x)) +
                       (static_cast<int64_t>(a.y - origin.y) * (a.y - origin.y));
    const int64_t db = (static_cast<int64_t>(b.x - origin.x) * (b.x - origin.x)) +
                       (static_cast<int64_t>(b.y - origin.y) * (b.y - origin.y));
    return da < db;
  });

  size_t m = 1;
  for (size_t i = 2; i < n; ++i) {
    while (m > 0 && Cross(pts[m - 1], pts[m], pts[i]) == 0) {
      --m;
    }
    ++m;
    pts[m] = pts[i];
  }
  pts.resize(m + 1);

  out.push_back(pts[0]);
  if (pts.size() <= 2) {
    if (pts.size() == 2) {
      out.push_back(pts[1]);
    }
    return;
  }
  out.push_back(pts[1]);
  for (size_t i = 2; i < pts.size(); ++i) {
    while (out.size() >= 2 && Cross(out[out.size() - 2], out.back(), pts[i]) <= 0) {
      out.pop_back();
    }
    out.push_back(pts[i]);
  }
}

void FloodFillComponent(const BinaryImage &img, int start_row, int start_col, std::vector<int> &label,
                        std::vector<Point> &component_pts) {
  const int rows = img.rows;
  const int cols = img.cols;
  component_pts.clear();
  std::vector<std::pair<int, int>> stack;
  stack.emplace_back(start_row, start_col);
  const size_t start_idx = img.Index(start_row, start_col);
  label[start_idx] = 1;
  while (!stack.empty()) {
    const auto [cur_r, cur_c] = stack.back();
    stack.pop_back();
    component_pts.push_back(Point{.x = cur_c, .y = cur_r});
    for (const auto &[dr, dc] : kFourNeighbors) {
      const int nr = cur_r + dr;
      const int nc = cur_c + dc;
      if (nr < 0 || nr >= rows || nc < 0 || nc >= cols) {
        continue;
      }
      const size_t nidx = img.Index(nr, nc);
      if (img.data[nidx] != 0 && label[nidx] == 0) {
        label[nidx] = 1;
        stack.emplace_back(nr, nc);
      }
    }
  }
}

std::vector<std::vector<Point>> ExtractComponents(const BinaryImage &img) {
  std::vector<std::vector<Point>> components;
  if (img.rows <= 0 || img.cols <= 0) {
    return components;
  }

  const size_t total = static_cast<size_t>(img.rows) * static_cast<size_t>(img.cols);
  std::vector<int> label(total, 0);
  std::vector<Point> component_pts;
  component_pts.reserve(total);

  for (int row = 0; row < img.rows; ++row) {
    for (int col = 0; col < img.cols; ++col) {
      const size_t idx = img.Index(row, col);
      if (img.data[idx] == 0 || label[idx] != 0) {
        continue;
      }
      FloodFillComponent(img, row, col, label, component_pts);
      components.push_back(component_pts);
    }
  }
  return components;
}

std::vector<int> FlattenComponents(const std::vector<std::vector<Point>> &comps) {
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

std::vector<std::vector<Point>> UnflattenComponents(const std::vector<int> &flat) {
  std::vector<std::vector<Point>> comps;
  std::size_t pos = 0;
  while (pos < flat.size()) {
    const int cnt = flat[pos++];
    std::vector<Point> comp;
    comp.reserve(static_cast<std::size_t>(std::max(0, cnt)));
    for (int i = 0; i < cnt && (pos + 1) < flat.size(); ++i) {
      comp.push_back(Point{.x = flat[pos], .y = flat[pos + 1]});
      pos += 2;
    }
    comps.push_back(std::move(comp));
  }
  return comps;
}

std::vector<int> MakeDispls(const std::vector<int> &counts) {
  std::vector<int> displs(counts.size(), 0);
  for (std::size_t i = 1; i < counts.size(); ++i) {
    displs[i] = displs[i - 1] + counts[i - 1];
  }
  return displs;
}

std::vector<int> MakeComponentCounts(int total_components, int proc_count) {
  std::vector<int> comp_counts(static_cast<std::size_t>(proc_count), 0);
  for (int i = 0; i < proc_count; ++i) {
    comp_counts[static_cast<std::size_t>(i)] =
        (total_components / proc_count) + ((i < (total_components % proc_count)) ? 1 : 0);
  }
  return comp_counts;
}

std::vector<int> FlattenDistributedComponents(const std::vector<std::vector<Point>> &all_components,
                                              const std::vector<int> &comp_counts, std::vector<int> &send_counts) {
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

    const std::vector<int> flat_part = FlattenComponents(part);
    send_counts[proc] = static_cast<int>(flat_part.size());
    flat_send.insert(flat_send.end(), flat_part.begin(), flat_part.end());
  }

  return flat_send;
}

HullList BuildLocalHulls(std::vector<std::vector<Point>> local_components) {
  HullList local_hulls(local_components.size());
  const int local_count = static_cast<int>(local_components.size());

#pragma omp parallel for default(none) shared(local_components, local_hulls, local_count) schedule(static)
  for (int i = 0; i < local_count; ++i) {
    const auto idx = static_cast<std::size_t>(i);
    std::vector<Point> pts = std::move(local_components[idx]);
    GrahamHull(pts, local_hulls[idx]);
  }

  return local_hulls;
}

std::vector<int> GatherFlatHulls(const std::vector<int> &flat_local_hulls, int rank, int size) {
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

HullList FlatToHullList(const std::vector<int> &flat) {
  return UnflattenComponents(flat);
}

HullList SolveAllMpi(const BinaryImage &img) {
  int rank = 0;
  int size = 1;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &size);

  std::vector<int> send_counts(static_cast<std::size_t>(size), 0);
  std::vector<int> send_displs(static_cast<std::size_t>(size), 0);
  std::vector<int> flat_send;

  if (rank == 0) {
    const auto all_components = ExtractComponents(img);
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
  const HullList local_hulls = BuildLocalHulls(std::move(local_components));
  const std::vector<int> flat_local_hulls = FlattenComponents(local_hulls);

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

  return FlatToHullList(flat_broadcast);
}

}  // namespace

KamalaginABinaryImageConvexHullALL::KamalaginABinaryImageConvexHullALL(const InType &in) {
  SetTypeOfTask(GetStaticTypeOfTask());
  GetInput() = in;
  GetOutput() = HullList{};
}

bool KamalaginABinaryImageConvexHullALL::ValidationImpl() {
  const auto &img = GetInput();
  if (img.rows < 0 || img.cols < 0) {
    return false;
  }
  if (img.rows == 0 || img.cols == 0) {
    return img.data.empty();
  }
  if (img.rows > 8192 || img.cols > 8192) {
    return false;
  }
  return (static_cast<size_t>(img.rows) * static_cast<size_t>(img.cols)) == img.data.size();
}

bool KamalaginABinaryImageConvexHullALL::PreProcessingImpl() {
  GetOutput().clear();
  return true;
}

bool KamalaginABinaryImageConvexHullALL::RunImpl() {
  GetOutput() = SolveAllMpi(GetInput());
  return true;
}

bool KamalaginABinaryImageConvexHullALL::PostProcessingImpl() {
  return true;
}

}  // namespace kamalagin_a_binary_image_convex_hull
