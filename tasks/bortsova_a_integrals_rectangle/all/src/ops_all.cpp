#include "bortsova_a_integrals_rectangle/all/include/ops_all.hpp"

#include <mpi.h>
#include <omp.h>

#include <algorithm>
#include <cstdint>
#include <vector>

#include "bortsova_a_integrals_rectangle/common/include/common.hpp"
#include "util/include/util.hpp"

namespace bortsova_a_integrals_rectangle {

BortsovaAIntegralsRectangleALL::BortsovaAIntegralsRectangleALL(const InType &in) {
  SetTypeOfTask(GetStaticTypeOfTask());
  GetInput() = in;
  GetOutput() = 0.0;
}

bool BortsovaAIntegralsRectangleALL::ValidationImpl() {
  const auto &input = GetInput();
  return input.func && !input.lower_bounds.empty() && input.lower_bounds.size() == input.upper_bounds.size() &&
         input.num_steps > 0;
}

bool BortsovaAIntegralsRectangleALL::PreProcessingImpl() {
  const auto &input = GetInput();
  func_ = input.func;
  num_steps_ = input.num_steps;
  dims_ = static_cast<int>(input.lower_bounds.size());

  midpoints_.resize(dims_);
  volume_ = 1.0;
  total_points_ = 1;

  for (int di = 0; di < dims_; di++) {
    double step = (input.upper_bounds[di] - input.lower_bounds[di]) / static_cast<double>(num_steps_);
    volume_ *= step;
    total_points_ *= num_steps_;

    midpoints_[di].resize(num_steps_);
    for (int si = 0; si < num_steps_; si++) {
      midpoints_[di][si] = input.lower_bounds[di] + ((si + 0.5) * step);
    }
  }

  return true;
}

double BortsovaAIntegralsRectangleALL::ComputePartialSum(int64_t begin, int64_t end) {
  std::vector<int> indices(dims_, 0);
  std::vector<double> point(dims_);

  int64_t temp = begin;
  for (int di = dims_ - 1; di >= 0; di--) {
    indices[di] = static_cast<int>(temp % num_steps_);
    temp /= num_steps_;
  }

  double local_sum = 0.0;
  for (int64_t pt = begin; pt < end; pt++) {
    for (int di = 0; di < dims_; di++) {
      point[di] = midpoints_[di][indices[di]];
    }
    local_sum += func_(point);

    for (int di = dims_ - 1; di >= 0; di--) {
      indices[di]++;
      if (indices[di] < num_steps_) {
        break;
      }
      indices[di] = 0;
    }
  }
  return local_sum;
}

bool BortsovaAIntegralsRectangleALL::RunImpl() {
  int rank = 0;
  int size = 1;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &size);

  int64_t chunk = total_points_ / size;
  int64_t remainder = total_points_ % size;
  int64_t mpi_begin = (rank * chunk) + std::min(static_cast<int64_t>(rank), remainder);
  int64_t mpi_end = mpi_begin + chunk + (static_cast<int64_t>(rank) < remainder ? 1 : 0);

  int num_threads = ppc::util::GetNumThreads();
  int64_t local_count = mpi_end - mpi_begin;
  std::vector<double> thread_sums(num_threads, 0.0);

#pragma omp parallel num_threads(num_threads) default(none) shared(thread_sums, mpi_begin, local_count, num_threads)
  {
    int tid = omp_get_thread_num();
    int64_t th_chunk = local_count / num_threads;
    int64_t th_rem = local_count % num_threads;
    int64_t th_begin = mpi_begin + (tid * th_chunk) + std::min(static_cast<int64_t>(tid), th_rem);
    int64_t th_end = th_begin + th_chunk + (static_cast<int64_t>(tid) < th_rem ? 1 : 0);
    thread_sums[tid] = ComputePartialSum(th_begin, th_end);
  }

  double local_sum = 0.0;
  for (int ti = 0; ti < num_threads; ti++) {
    local_sum += thread_sums[ti];
  }

  double global_sum = 0.0;
  MPI_Reduce(&local_sum, &global_sum, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
  MPI_Bcast(&global_sum, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);

  GetOutput() = global_sum * volume_;
  return true;
}

bool BortsovaAIntegralsRectangleALL::PostProcessingImpl() {
  return true;
}

}  // namespace bortsova_a_integrals_rectangle
