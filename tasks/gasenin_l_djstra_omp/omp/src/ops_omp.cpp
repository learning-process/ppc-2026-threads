#include "gasenin_l_djstra_omp/omp/include/ops_omp.hpp"

#include <omp.h>

#include <algorithm>
#include <cstdint>
#include <cstdlib>
#include <limits>
#include <vector>

#include "gasenin_l_djstra_omp/common/include/common.hpp"

namespace gasenin_l_djstra_omp {

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

  auto &dist = dist_;
  auto &visited = visited_;

  int num_threads = 1;

#pragma omp parallel default(none) shared(num_threads)
  {
#pragma omp single
    num_threads = omp_get_num_threads();
  }

  std::vector<InType> local_min(num_threads, inf);
  std::vector<InType> local_vertex(num_threads, -1);

  InType global_vertex = -1;

#pragma omp parallel default(none) shared(n, inf, dist, visited, local_min, local_vertex, global_vertex, num_threads)
  {
    const int thread_id = omp_get_thread_num();

    for (int iteration = 0; iteration < n; ++iteration) {
      InType thread_min = inf;
      InType thread_vertex = -1;

#pragma omp for nowait
      for (int index = 0; index < n; ++index) {
        if (visited_[index] == 0 && dist_[index] < thread_min) {
          thread_min = dist_[index];
          thread_vertex = index;
        }
      }

      local_min[thread_id] = thread_min;
      local_vertex[thread_id] = thread_vertex;

#pragma omp barrier

#pragma omp single
      {
        global_vertex = -1;
        InType global_min = inf;

        for (int thread_idx = 0; thread_idx < num_threads; ++thread_idx) {
          if (local_min[thread_idx] < global_min) {
            global_min = local_min[thread_idx];
            global_vertex = local_vertex[thread_idx];
          }
        }

        if (global_vertex != -1 && global_min != inf) {
          visited_[global_vertex] = 1;
        }
      }

#pragma omp barrier

      if (global_vertex == -1) {
        break;
      }

#pragma omp for
      for (int vertex = 0; vertex < n; ++vertex) {
        if (visited_[vertex] == 0 && vertex != global_vertex) {
          const InType weight = std::abs(global_vertex - vertex);

          if (dist_[global_vertex] != inf) {
            const InType new_dist = dist_[global_vertex] + weight;
            dist_[vertex] = std::min(dist_[vertex], new_dist);
          }
        }
      }
    }
  }

  int64_t total_sum = 0;

#pragma omp parallel for reduction(+ : total_sum) default(none) shared(n, dist, inf)
  for (int index = 0; index < n; ++index) {
    if (dist_[index] != inf) {
      total_sum += dist_[index];
    }
  }

  GetOutput() = static_cast<OutType>(total_sum);
  return true;
}

bool GaseninLDjstraOMP::PostProcessingImpl() {
  return true;
}

}  // namespace gasenin_l_djstra_omp
