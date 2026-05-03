#include "gasenin_l_djstra/omp/include/ops_omp.hpp"

#include <omp.h>

#include <algorithm>
#include <cstdint>
#include <cstdlib>
#include <limits>
#include <vector>

#include "gasenin_l_djstra/common/include/common.hpp"

namespace gasenin_l_djstra {

namespace {

int GetNumThreads() {
  int num_threads = 1;

#pragma omp parallel default(none) shared(num_threads)
  {
#pragma omp single
    num_threads = omp_get_num_threads();
  }

  return num_threads;
}

InType FindGlobalVertexOMP(const InType n, const InType inf, std::vector<InType> &dist, std::vector<char> &visited,
                           const int num_threads) {
  std::vector<InType> local_min(num_threads, inf);
  std::vector<InType> local_vertex(num_threads, -1);

  InType global_vertex = -1;

#pragma omp parallel default(none) shared(n, inf, dist, visited, local_min, local_vertex, global_vertex, num_threads)
  {
    const int thread_id = omp_get_thread_num();

    InType thread_min = inf;
    InType thread_vertex = -1;

#pragma omp for nowait
    for (int index = 0; index < n; ++index) {
      if (visited[index] == 0 && dist[index] < thread_min) {
        thread_min = dist[index];
        thread_vertex = index;
      }
    }

    local_min[thread_id] = thread_min;
    local_vertex[thread_id] = thread_vertex;

#pragma omp barrier

#pragma omp single
    {
      InType global_min = inf;

      for (int thread_idx = 0; thread_idx < num_threads; ++thread_idx) {
        if (local_min[thread_idx] < global_min) {
          global_min = local_min[thread_idx];
          global_vertex = local_vertex[thread_idx];
        }
      }

      if (global_vertex != -1 && global_min != inf) {
        visited[global_vertex] = 1;
      } else {
        global_vertex = -1;
      }
    }
  }

  return global_vertex;
}

void RelaxEdgesOMP(const InType n, const InType inf, const InType global_vertex, std::vector<InType> &dist,
                   std::vector<char> &visited) {
#pragma omp parallel for default(none) shared(n, inf, global_vertex, dist, visited)
  for (int vertex = 0; vertex < n; ++vertex) {
    if (visited[vertex] == 0 && vertex != global_vertex) {
      const InType weight = std::abs(global_vertex - vertex);

      if (dist[global_vertex] != inf) {
        const InType new_dist = dist[global_vertex] + weight;
        dist[vertex] = std::min(dist[vertex], new_dist);
      }
    }
  }
}

int64_t CalculateTotalSumOMP(const InType n, const InType inf, std::vector<InType> &dist) {
  int64_t total_sum = 0;

#pragma omp parallel for reduction(+ : total_sum) default(none) shared(n, inf, dist)
  for (int i = 0; i < n; ++i) {
    if (dist[i] != inf) {
      total_sum += dist[i];
    }
  }

  return total_sum;
}

}  // namespace

GaseninLDjstraOMP::GaseninLDjstraOMP(const InType &in) {
  SetTypeOfTask(GetStaticTypeOfTask());
  GetInput() = in;
  GetOutput() = 0;
}

bool GaseninLDjstraOMP::ValidationImpl() {
  return GetInput() > 0;
}

bool GaseninLDjstraOMP::PreProcessingImpl() {
  const InType n = GetInput();
  const InType inf = std::numeric_limits<InType>::max();

  dist_.assign(n, inf);
  visited_.assign(n, 0);

  dist_[0] = 0;
  return true;
}

bool GaseninLDjstraOMP::RunImpl() {
  const InType n = GetInput();
  const InType inf = std::numeric_limits<InType>::max();

  const int num_threads = GetNumThreads();

  for (int iteration = 0; iteration < n; ++iteration) {
    InType global_vertex = FindGlobalVertexOMP(n, inf, dist_, visited_, num_threads);

    if (global_vertex == -1) {
      break;
    }

    RelaxEdgesOMP(n, inf, global_vertex, dist_, visited_);
  }

  const int64_t total_sum = CalculateTotalSumOMP(n, inf, dist_);

  GetOutput() = static_cast<OutType>(total_sum);
  return true;
}

bool GaseninLDjstraOMP::PostProcessingImpl() {
  return true;
}

}  // namespace gasenin_l_djstra
