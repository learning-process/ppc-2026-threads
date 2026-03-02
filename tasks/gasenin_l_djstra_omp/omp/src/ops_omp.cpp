#include "gasenin_l_djstra_omp/omp/include/ops_omp.hpp"

#include <omp.h>

#include <algorithm>
#include <limits>
#include <vector>

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
  const InType INF = std::numeric_limits<InType>::max();

  dist_.assign(n, INF);
  visited_.assign(n, 0);

  dist_[0] = 0;
  return true;
}

bool GaseninLDjstraOMP::RunImpl() {
  const InType n = GetInput();
  const InType INF = std::numeric_limits<InType>::max();

  int num_threads = 1;
#pragma omp parallel
  {
#pragma omp single
    num_threads = omp_get_num_threads();
  }

  std::vector<InType> local_min(num_threads, INF);
  std::vector<InType> local_u(num_threads, -1);

  InType global_u = -1;
  bool done = false;

#pragma omp parallel
  {
    int thread_id = omp_get_thread_num();

    for (int count = 0; count < n && !done; ++count) {
      InType my_min = INF;
      InType my_u = -1;

#pragma omp for nowait
      for (int i = 0; i < n; ++i) {
        if (!visited_[i] && dist_[i] < my_min) {
          my_min = dist_[i];
          my_u = i;
        }
      }

      local_min[thread_id] = my_min;
      local_u[thread_id] = my_u;

#pragma omp barrier

#pragma omp single
      {
        global_u = -1;
        InType global_min = INF;
        for (int t = 0; t < num_threads; ++t) {
          if (local_min[t] < global_min) {
            global_min = local_min[t];
            global_u = local_u[t];
          }
        }
        if (global_u == -1 || global_min == INF) {
          done = true;
        } else {
          visited_[global_u] = 1;
        }
      }

      if (done) {
        break;
      }

#pragma omp for
      for (int v = 0; v < n; ++v) {
        if (!visited_[v] && v != global_u) {
          InType weight = std::abs(global_u - v);

          if (dist_[global_u] != INF) {
            InType new_dist = dist_[global_u] + weight;
            if (new_dist < dist_[v]) {
              dist_[v] = new_dist;
            }
          }
        }
      }
    }
  }

  long long total_sum = 0;
#pragma omp parallel for reduction(+ : total_sum)
  for (int i = 0; i < n; ++i) {
    if (dist_[i] != INF) {
      total_sum += dist_[i];
    }
  }
  GetOutput() = static_cast<OutType>(total_sum);

  return true;
}

bool GaseninLDjstraOMP::PostProcessingImpl() {
  return true;
}

}  // namespace gasenin_l_djstra_omp
