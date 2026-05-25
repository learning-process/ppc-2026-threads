#include "egorova_l_binary_convex_hull/all/include/ops_all.hpp"

#include <mpi.h>
#include <omp.h>

#include <algorithm>
#include <cstddef>
#include <cstdint>
#include <functional>
#include <limits>
#include <queue>
#include <utility>
#include <vector>

#include "egorova_l_binary_convex_hull/common/include/common.hpp"

namespace egorova_l_binary_convex_hull {

namespace {

struct Pix {
  int px, py;
};

// Сериализация оболочек для MPI
std::vector<int> SerializeHulls(const std::vector<std::vector<Point>> &hulls) {
  std::vector<int> buf;
  buf.push_back(static_cast<int>(hulls.size()));
  for (const auto &hull : hulls) {
    buf.push_back(static_cast<int>(hull.size()));
    for (const auto &pp : hull) {
      buf.push_back(pp.x);
      buf.push_back(pp.y);
    }
  }
  return buf;
}

std::vector<std::vector<Point>> DeserializeHulls(const std::vector<int> &buf) {
  std::vector<std::vector<Point>> hulls;
  if (buf.empty()) {
    return hulls;
  }
  size_t pos = 0;
  const int nn = buf[pos++];
  hulls.resize(nn);
  for (int ii = 0; ii < nn; ++ii) {
    const int sz = buf[pos++];
    hulls[ii].reserve(sz);
    for (int jj = 0; jj < sz; ++jj) {
      hulls[ii].push_back({buf[pos], buf[pos + 1]});
      pos += 2;
    }
  }
  return hulls;
}

// Проверка соседа в BFS (локальная полоса)
inline void TryPushLocal(const std::vector<uint8_t> &image, int width, int local_rows, int row_start, int nx, int ny,
                         std::vector<uint8_t> &visited, std::queue<Pix> &bfs) {
  if (nx < 0 || nx >= width || ny < 0 || ny >= local_rows) {
    return;
  }
  const size_t local_ni = ((static_cast<size_t>(ny) * static_cast<size_t>(width)) + static_cast<size_t>(nx));
  const size_t global_ni =
      ((static_cast<size_t>(ny + row_start) * static_cast<size_t>(width)) + static_cast<size_t>(nx));
  if (image[global_ni] != 0 && visited[local_ni] == 0) {
    visited[local_ni] = 1;
    bfs.push({nx, ny});
  }
}

// Один BFS, возвращает компоненту (в глобальных координатах)
std::vector<Point> BfsOneComponent(const std::vector<uint8_t> &image, int width, int local_rows, int row_start,
                                   int start_x, int start_y, std::vector<uint8_t> &visited) {
  std::vector<Point> comp;
  comp.reserve(256);
  std::queue<Pix> bfs;
  bfs.push({start_x, start_y});
  while (!bfs.empty()) {
    const auto cur = bfs.front();
    bfs.pop();
    comp.push_back({cur.px, cur.py + row_start});
    TryPushLocal(image, width, local_rows, row_start, cur.px + 1, cur.py, visited, bfs);
    TryPushLocal(image, width, local_rows, row_start, cur.px - 1, cur.py, visited, bfs);
    TryPushLocal(image, width, local_rows, row_start, cur.px, cur.py + 1, visited, bfs);
    TryPushLocal(image, width, local_rows, row_start, cur.px, cur.py - 1, visited, bfs);
  }
  return comp;
}

std::vector<std::vector<Point>> FindLocalComponents(const std::vector<uint8_t> &image, int width, int local_rows,
                                                    int row_start) {
  std::vector<std::vector<Point>> comps;
  std::vector<uint8_t> visited(((static_cast<size_t>(local_rows) * static_cast<size_t>(width))), 0);

  for (int ly = 0; ly < local_rows; ++ly) {
    const int gy = row_start + ly;
    for (int col = 0; col < width; ++col) {
      const size_t local_idx = ((static_cast<size_t>(ly) * static_cast<size_t>(width)) + static_cast<size_t>(col));
      const size_t global_idx = ((static_cast<size_t>(gy) * static_cast<size_t>(width)) + static_cast<size_t>(col));
      if (image[global_idx] == 0 || visited[local_idx] != 0) {
        continue;
      }
      visited[local_idx] = 1;
      comps.push_back(BfsOneComponent(image, width, local_rows, row_start, col, ly, visited));
    }
  }
  return comps;
}

// Вычисление bounding box'ов всех оболочек
struct BBox {
  std::vector<int> min_x, max_x, min_y, max_y;
};

BBox ComputeBBoxes(const std::vector<std::vector<Point>> &hulls) {
  const int nn = static_cast<int>(hulls.size());
  BBox bb;
  bb.min_x.assign(nn, std::numeric_limits<int>::max());
  bb.max_x.assign(nn, std::numeric_limits<int>::min());
  bb.min_y.assign(nn, std::numeric_limits<int>::max());
  bb.max_y.assign(nn, std::numeric_limits<int>::min());
  for (int ii = 0; ii < nn; ++ii) {
    for (const auto &pp : hulls[ii]) {
      bb.min_x[ii] = std::min(bb.min_x[ii], pp.x);
      bb.max_x[ii] = std::max(bb.max_x[ii], pp.x);
      bb.min_y[ii] = std::min(bb.min_y[ii], pp.y);
      bb.max_y[ii] = std::max(bb.max_y[ii], pp.y);
    }
  }
  return bb;
}

// Union-Find для объединения смежных по Y и пересекающихся по X оболочек
void UnionAdjacent(int nn, const BBox &bb, std::vector<int> &parent, const std::function<int(int)> &find) {
  for (int ii = 0; ii < nn; ++ii) {
    for (int jj = ii + 1; jj < nn; ++jj) {
      if (find(ii) == find(jj)) {
        continue;
      }
      const bool y_adj = (bb.max_y[ii] + 1 == bb.min_y[jj]) || (bb.max_y[jj] + 1 == bb.min_y[ii]);
      const bool x_ovl = (bb.min_x[ii] <= bb.max_x[jj]) && (bb.min_x[jj] <= bb.max_x[ii]);
      if (y_adj && x_ovl) {
        parent[find(ii)] = find(jj);
      }
    }
  }
}

std::vector<std::vector<Point>> MergeHulls(const std::vector<std::vector<Point>> &all_hulls) {
  const int nn = static_cast<int>(all_hulls.size());
  std::vector<int> parent(nn);
  for (int ii = 0; ii < nn; ++ii) {
    parent[ii] = ii;
  }
  std::function<int(int)> find = [&](int ii) -> int { return parent[ii] == ii ? ii : (parent[ii] = find(parent[ii])); };

  const BBox bb = ComputeBBoxes(all_hulls);
  UnionAdjacent(nn, bb, parent, find);

  std::vector<std::vector<Point>> groups(nn);
  for (int ii = 0; ii < nn; ++ii) {
    const int root = find(ii);
    for (const auto &pp : all_hulls[ii]) {
      groups[root].push_back(pp);
    }
  }
  return groups;
}

// Сборка всех оболочек со всех процессов через Allgatherv
std::vector<std::vector<Point>> GatherAllHulls(const std::vector<std::vector<Point>> &local_hulls, int p_count) {
  std::vector<int> s_buf = SerializeHulls(local_hulls);
  const int s_size = static_cast<int>(s_buf.size());

  std::vector<int> r_sizes(p_count, 0);
  MPI_Allgather(&s_size, 1, MPI_INT, r_sizes.data(), 1, MPI_INT, MPI_COMM_WORLD);

  std::vector<int> r_displs(p_count, 0);
  for (int ii = 1; ii < p_count; ++ii) {
    r_displs[ii] = r_displs[ii - 1] + r_sizes[ii - 1];
  }
  const int total_size = r_displs[p_count - 1] + r_sizes[p_count - 1];

  std::vector<int> r_buf(total_size);
  MPI_Allgatherv(s_buf.data(), s_size, MPI_INT, r_buf.data(), r_sizes.data(), r_displs.data(), MPI_INT, MPI_COMM_WORLD);

  std::vector<std::vector<Point>> all_hulls;
  all_hulls.reserve(64);
  for (int ii = 0; ii < p_count; ++ii) {
    if (r_sizes[ii] == 0) {
      continue;
    }
    std::vector<int> ph(r_buf.begin() + r_displs[ii], r_buf.begin() + r_displs[ii] + r_sizes[ii]);
    auto dec = DeserializeHulls(ph);
    for (auto &hh : dec) {
      all_hulls.push_back(std::move(hh));
    }
  }
  return all_hulls;
}

// BFS для FindComponents (без полос)
inline void TryPushGlobal(const std::vector<uint8_t> &image, int width, int height, int nx, int ny,
                          std::vector<uint8_t> &visited, std::queue<Pix> &bfs) {
  if (nx < 0 || nx >= width || ny < 0 || ny >= height) {
    return;
  }
  const size_t ni = ((static_cast<size_t>(ny) * static_cast<size_t>(width)) + static_cast<size_t>(nx));
  if (image[ni] != 0 && visited[ni] == 0) {
    visited[ni] = 1;
    bfs.push({nx, ny});
  }
}

std::vector<Point> BfsComponentGlobal(const std::vector<uint8_t> &image, int width, int height, int start_x,
                                      int start_y, std::vector<uint8_t> &visited) {
  std::vector<Point> comp;
  std::queue<Pix> bfs;
  bfs.push({start_x, start_y});
  while (!bfs.empty()) {
    const auto cur = bfs.front();
    bfs.pop();
    comp.emplace_back(cur.px, cur.py);
    TryPushGlobal(image, width, height, cur.px + 1, cur.py, visited, bfs);
    TryPushGlobal(image, width, height, cur.px - 1, cur.py, visited, bfs);
    TryPushGlobal(image, width, height, cur.px, cur.py + 1, visited, bfs);
    TryPushGlobal(image, width, height, cur.px, cur.py - 1, visited, bfs);
  }
  return comp;
}

// BFS с labels (для ProcessComponent)
inline void TryPushLabels(const std::vector<uint8_t> &image, int width, int height, int nx, int ny, int label,
                          std::vector<int> &labels, std::queue<Pix> &bfs) {
  if (nx < 0 || nx >= width || ny < 0 || ny >= height) {
    return;
  }
  const size_t ni = ((static_cast<size_t>(ny) * static_cast<size_t>(width)) + static_cast<size_t>(nx));
  if (image[ni] != 0 && labels[ni] == 0) {
    labels[ni] = label;
    bfs.push({nx, ny});
  }
}

}  // namespace

BinaryConvexHullALL::BinaryConvexHullALL(const InType &in) {
  SetTypeOfTask(GetStaticTypeOfTask());
  GetInput() = in;
}

bool BinaryConvexHullALL::ValidationImpl() {
  int rank = 0;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  if (rank == 0) {
    const auto &in = GetInput();
    return in.width > 0 && in.height > 0 && !in.data.empty();
  }
  return true;
}

bool BinaryConvexHullALL::PreProcessingImpl() {
  GetOutput().clear();
  return true;
}

bool BinaryConvexHullALL::RunImpl() {
  int rank = 0;
  int p_count = 1;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &p_count);

  // Bcast размеров и изображения
  int width = 0;
  int height = 0;
  if (rank == 0) {
    width = GetInput().width;
    height = GetInput().height;
  }
  MPI_Bcast(&width, 1, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Bcast(&height, 1, MPI_INT, 0, MPI_COMM_WORLD);

  if (width <= 0 || height <= 0) {
    GetOutput().clear();
    return true;
  }

  std::vector<uint8_t> image;
  if (rank == 0) {
    image = GetInput().data;
  } else {
    image.resize((static_cast<size_t>(width) * static_cast<size_t>(height)));
  }
  MPI_Bcast(image.data(), static_cast<int>(image.size()), MPI_BYTE, 0, MPI_COMM_WORLD);

  // Распределение строк
  const int chunk = height / p_count;
  const int rem = height % p_count;
  const int row_start = ((rank * chunk) + std::min(rank, rem));
  const int local_rows = chunk + ((rank < rem) ? 1 : 0);

  // BFS локально + OpenMP для оболочек
  std::vector<std::vector<Point>> local_comps = FindLocalComponents(image, width, local_rows, row_start);
  std::vector<std::vector<Point>> local_hulls(local_comps.size());
  const int n_local = static_cast<int>(local_comps.size());
#pragma omp parallel for schedule(dynamic, 1) default(none) shared(local_comps, local_hulls, n_local)
  for (int ii = 0; ii < n_local; ++ii) {
    BuildConvexHull(local_comps[ii], local_hulls[ii]);
  }

  // Allgatherv → каждый процесс получает все оболочки
  std::vector<std::vector<Point>> all_hulls = GatherAllHulls(local_hulls, p_count);
  if (all_hulls.empty()) {
    GetOutput().clear();
    return true;
  }

  // Union-find + финальные оболочки
  std::vector<std::vector<Point>> groups = MergeHulls(all_hulls);
  GetOutput().clear();
  for (auto &group : groups) {
    if (!group.empty()) {
      std::vector<Point> hh;
      BuildConvexHull(group, hh);
      GetOutput().push_back(std::move(hh));
    }
  }
  return true;
}

bool BinaryConvexHullALL::PostProcessingImpl() {
  return true;
}

int64_t BinaryConvexHullALL::CrossProduct(const Point &a, const Point &b, const Point &c) {
  return (((static_cast<int64_t>(b.x - a.x) * static_cast<int64_t>(c.y - a.y)) -
           (static_cast<int64_t>(b.y - a.y) * static_cast<int64_t>(c.x - a.x))));
}

void BinaryConvexHullALL::BuildConvexHull(std::vector<Point> &points, std::vector<Point> &hull) {
  hull.clear();
  if (points.empty()) {
    return;
  }
  std::ranges::sort(points, [](const Point &a, const Point &b) { return a.x < b.x || (a.x == b.x && a.y < b.y); });
  const auto last_unique = std::ranges::unique(points).begin();
  points.erase(last_unique, points.end());
  if (points.size() < 3) {
    hull.clear();
    for (const auto &pp : points) {
      hull.push_back(pp);
    }
    return;
  }
  std::vector<Point> res;
  res.reserve(points.size() + 1);
  for (const auto &pp : points) {
    while (res.size() >= 2 && CrossProduct(res[res.size() - 2], res.back(), pp) <= 0) {
      res.pop_back();
    }
    res.push_back(pp);
  }
  const size_t lower_size = res.size();
  for (int ii = static_cast<int>(points.size()) - 2; ii >= 0; --ii) {
    while (res.size() > lower_size && CrossProduct(res[res.size() - 2], res.back(), points[ii]) <= 0) {
      res.pop_back();
    }
    res.push_back(points[ii]);
  }
  if (res.size() > 1) {
    res.pop_back();
  }
  hull.swap(res);
}

std::vector<std::vector<Point>> BinaryConvexHullALL::FindComponents(const std::vector<uint8_t> &image, int width,
                                                                    int height) {
  std::vector<std::vector<Point>> components;
  std::vector<uint8_t> visited(((static_cast<size_t>(width) * static_cast<size_t>(height))), 0);
  for (int row = 0; row < height; ++row) {
    for (int col = 0; col < width; ++col) {
      const size_t idx = ((static_cast<size_t>(row) * static_cast<size_t>(width)) + static_cast<size_t>(col));
      if (image[idx] == 0 || visited[idx] != 0) {
        continue;
      }
      visited[idx] = 1;
      components.push_back(BfsComponentGlobal(image, width, height, col, row, visited));
    }
  }
  return components;
}

void BinaryConvexHullALL::ProcessComponent(const std::vector<uint8_t> &image, int width, int height, int start_x,
                                           int start_y, int label, std::vector<int> &labels,
                                           std::vector<Point> &component_points) {
  if (start_x < 0 || start_x >= width || start_y < 0 || start_y >= height) {
    component_points.clear();
    return;
  }
  const size_t si = ((static_cast<size_t>(start_y) * static_cast<size_t>(width)) + static_cast<size_t>(start_x));
  labels[si] = label;
  std::queue<Pix> bfs;
  bfs.push({start_x, start_y});
  component_points.clear();
  component_points.reserve(1000);
  while (!bfs.empty()) {
    const auto cur = bfs.front();
    bfs.pop();
    component_points.emplace_back(cur.px, cur.py);
    TryPushLabels(image, width, height, cur.px + 1, cur.py, label, labels, bfs);
    TryPushLabels(image, width, height, cur.px - 1, cur.py, label, labels, bfs);
    TryPushLabels(image, width, height, cur.px, cur.py + 1, label, labels, bfs);
    TryPushLabels(image, width, height, cur.px, cur.py - 1, label, labels, bfs);
  }
}

}  // namespace egorova_l_binary_convex_hull
