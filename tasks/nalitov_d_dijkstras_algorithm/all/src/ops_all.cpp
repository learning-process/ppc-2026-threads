#include "nalitov_d_dijkstras_algorithm/all/include/ops_all.hpp"

#include <mpi.h>
#include <omp.h>

#include <algorithm>
#include <array>
#include <atomic>
#include <cstddef>
#include <cstdint>
#include <utility>
#include <vector>

#include "nalitov_d_dijkstras_algorithm/common/include/common.hpp"

namespace nalitov_d_dijkstras_algorithm {

namespace {

struct NodePayload {
  int cost;
  int vtx;
};

NodePayload FindLocalMin(const std::vector<Cost> &dist, const std::vector<char> &visited, int start_idx, int count,
                         int threads) {
  int best_c = kInf;
  int best_v = -1;

#pragma omp parallel num_threads(threads) default(none) shared(dist, visited, start_idx, count, best_c, best_v)
  {
    int thr_c = kInf;
    int thr_v = -1;

#pragma omp for nowait
    for (int i = 0; i < count; ++i) {
      if (visited[static_cast<size_t>(i)] == 0 && dist[static_cast<size_t>(i)] < thr_c) {
        thr_c = dist[static_cast<size_t>(i)];
        thr_v = start_idx + i;
      }
    }

#pragma omp critical
    {
      if (thr_c < best_c || (thr_c == best_c && thr_v != -1 && (best_v == -1 || thr_v < best_v))) {
        best_c = thr_c;
        best_v = thr_v;
      }
    }
  }

  return {.cost = best_c, .vtx = best_v};
}

void UpdateNeighborDistance(int target, int weight, int global_cost, int l_start, int l_count,
                            std::vector<Cost> &local_dist, const std::vector<char> &local_visited) {
  if (target >= l_start && target < l_start + l_count) {
    int local_idx = target - l_start;
    if (local_visited[static_cast<size_t>(local_idx)] == 0) {
      int new_dist = global_cost + weight;
      std::atomic_ref<int> target_ref(local_dist[static_cast<size_t>(local_idx)]);

      int old_dist = target_ref.load(std::memory_order_relaxed);
      while (new_dist < old_dist) {
        if (target_ref.compare_exchange_weak(old_dist, new_dist, std::memory_order_relaxed)) {
          break;
        }
      }
    }
  }
}

void RelaxNeighbors(const std::vector<std::pair<int, int>> &neighbors, const NodePayload &global, int l_start,
                    int l_count, std::vector<Cost> &local_dist, const std::vector<char> &local_visited, int threads) {
  const size_t sz = neighbors.size();
#pragma omp parallel for num_threads(threads) default(none) \
    shared(neighbors, sz, global, l_start, l_count, local_dist, local_visited)
  for (size_t i = 0; i < sz; ++i) {
    int target = neighbors[i].first;
    int weight = neighbors[i].second;
    UpdateNeighborDistance(target, weight, global.cost, l_start, l_count, local_dist, local_visited);
  }
}

int64_t CalculateLocalSum(const std::vector<Cost> &local_dist, int l_count, int threads) {
  int64_t local_sum = 0;
#pragma omp parallel for num_threads(threads) default(none) reduction(+ : local_sum) shared(local_dist, l_count)
  for (int i = 0; i < l_count; ++i) {
    if (local_dist[static_cast<size_t>(i)] != kInf) {
      local_sum += local_dist[static_cast<size_t>(i)];
    }
  }
  return local_sum;
}

}  // namespace

NalitovDDijkstrasAlgorithmALL::NalitovDDijkstrasAlgorithmALL(const InType &in) {
  SetTypeOfTask(GetStaticTypeOfTask());
  GetInput() = in;
  GetOutput() = 0;
}

bool NalitovDDijkstrasAlgorithmALL::ValidationImpl() {
  MPI_Comm_rank(MPI_COMM_WORLD, &rank_);
  int is_valid = 1;

  if (rank_ == 0) {
    const auto &in = GetInput();
    if (in.n <= 0 || in.n > 10000 || in.source < 0 || in.source >= in.n) {
      is_valid = 0;
    } else {
      for (const auto &a : in.arcs) {
        if (a.from < 0 || a.to < 0 || a.from >= in.n || a.to >= in.n || a.weight < 0) {
          is_valid = 0;
          break;
        }
      }
    }
  }

  MPI_Bcast(&is_valid, 1, MPI_INT, 0, MPI_COMM_WORLD);
  return is_valid == 1;
}

bool NalitovDDijkstrasAlgorithmALL::PreProcessingImpl() {
  MPI_Comm_rank(MPI_COMM_WORLD, &rank_);
  MPI_Comm_size(MPI_COMM_WORLD, &size_);

  std::array<int, 3> meta{};
  if (rank_ == 0) {
    meta[0] = GetInput().n;
    meta[1] = GetInput().source;
    meta[2] = static_cast<int>(GetInput().arcs.size());
  }
  MPI_Bcast(meta.data(), 3, MPI_INT, 0, MPI_COMM_WORLD);

  n_ = meta[0];
  source_ = meta[1];
  int arc_count = meta[2];

  std::vector<int> arc_buf(static_cast<size_t>(arc_count) * 3);
  if (rank_ == 0) {
    for (int i = 0; i < arc_count; ++i) {
      size_t idx = static_cast<size_t>(i) * 3;
      arc_buf[idx + 0] = GetInput().arcs[i].from;
      arc_buf[idx + 1] = GetInput().arcs[i].to;
      arc_buf[idx + 2] = GetInput().arcs[i].weight;
    }
  }
  if (arc_count > 0) {
    MPI_Bcast(arc_buf.data(), arc_count * 3, MPI_INT, 0, MPI_COMM_WORLD);
  }

  graph_.assign(n_, {});
  for (int i = 0; i < arc_count; ++i) {
    size_t idx = static_cast<size_t>(i) * 3;
    graph_[static_cast<size_t>(arc_buf[idx])].emplace_back(arc_buf[idx + 1], arc_buf[idx + 2]);
  }

  const int base = n_ / size_;
  const int rem = n_ % size_;
  local_start_ = (rank_ * base) + std::min(rank_, rem);
  local_count_ = base + (rank_ < rem ? 1 : 0);

  local_dist_.assign(local_count_, kInf);
  local_visited_.assign(local_count_, 0);

  if (source_ >= local_start_ && source_ < local_start_ + local_count_) {
    local_dist_[source_ - local_start_] = 0;
  }

  return true;
}

bool NalitovDDijkstrasAlgorithmALL::RunImpl() {
  const int threads = omp_get_max_threads();
  const int l_start = local_start_;
  const int l_count = local_count_;
  const int n_total = n_;

  bool all_done = false;

  for (int step = 0; step < n_total; ++step) {
    NodePayload local = FindLocalMin(local_dist_, local_visited_, l_start, l_count, threads);

    NodePayload global{};
    MPI_Allreduce(&local, &global, 1, MPI_2INT, MPI_MINLOC, MPI_COMM_WORLD);

    if (global.vtx == -1 || global.cost == kInf) {
      all_done = true;
    }

    if (!all_done) {
      if (global.vtx >= l_start && global.vtx < l_start + l_count) {
        local_visited_[static_cast<size_t>(global.vtx - l_start)] = 1;
      }

      RelaxNeighbors(graph_[static_cast<size_t>(global.vtx)], global, l_start, l_count, local_dist_, local_visited_,
                     threads);
    }
  }

  int64_t local_sum = CalculateLocalSum(local_dist_, l_count, threads);
  int64_t global_sum = 0;
  MPI_Reduce(&local_sum, &global_sum, 1, MPI_INT64_T, MPI_SUM, 0, MPI_COMM_WORLD);

  if (rank_ == 0) {
    GetOutput() = global_sum;
  }

  return true;
}

bool NalitovDDijkstrasAlgorithmALL::PostProcessingImpl() {
  OutType result = GetOutput();
  MPI_Bcast(&result, 1, MPI_INT64_T, 0, MPI_COMM_WORLD);
  if (rank_ != 0) {
    GetOutput() = result;
  }
  return true;
}

}  // namespace nalitov_d_dijkstras_algorithm
