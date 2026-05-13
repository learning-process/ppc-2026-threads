#include "telnov_a_integral_rectangle/all/include/ops_all.hpp"

#include <mpi.h>
#include <omp.h>

#include <algorithm>
#include <cstdint>
#include <functional>
#include <thread>
#include <vector>

#include "oneapi/tbb/blocked_range.h"
#include "oneapi/tbb/parallel_reduce.h"
#include "telnov_a_integral_rectangle/common/include/common.hpp"
#include "util/include/util.hpp"

namespace telnov_a_integral_rectangle {

namespace {

int64_t PowInt(int base, int degree) {
  int64_t result = 1;
  for (int i = 0; i < degree; ++i) {
    result *= base;
  }
  return result;
}

double CalculatePointValue(int64_t index, int n, int dimensions, double h) {
  double value = 0.0;

  for (int dim = 0; dim < dimensions; ++dim) {
    const int coordinate_index = static_cast<int>(index % n);
    index /= n;

    const double x = ((static_cast<double>(coordinate_index) + 0.5) * h);
    value += x;
  }

  return value;
}

double CalculateRange(int64_t begin, int64_t end, int n, int dimensions, double h) {
  return oneapi::tbb::parallel_reduce(
      oneapi::tbb::blocked_range<int64_t>(begin, end), 0.0,
      [n, dimensions, h](const oneapi::tbb::blocked_range<int64_t> &range, double local_sum) {
    for (int64_t index = range.begin(); index != range.end(); ++index) {
      local_sum += CalculatePointValue(index, n, dimensions, h);
    }
    return local_sum;
  }, std::plus<>());
}

}  // namespace

TelnovAIntegralRectangleALL::TelnovAIntegralRectangleALL(const InType &in) {
  SetTypeOfTask(GetStaticTypeOfTask());
  GetInput() = in;
  GetOutput() = 0.0;
}

bool TelnovAIntegralRectangleALL::ValidationImpl() {
  return GetInput().first > 0 && GetInput().second > 0;
}

bool TelnovAIntegralRectangleALL::PreProcessingImpl() {
  GetOutput() = 0.0;
  return true;
}

bool TelnovAIntegralRectangleALL::RunImpl() {
  const int n = GetInput().first;
  const int dimensions = GetInput().second;

  const double h = 1.0 / static_cast<double>(n);
  const int64_t total_points = PowInt(n, dimensions);

  int rank = 0;
  int size = 1;

  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &size);

  const int64_t base_block = total_points / size;
  const int64_t remainder = total_points % size;

  const int64_t rank_begin = (static_cast<int64_t>(rank) * base_block) + std::min<int64_t>(rank, remainder);
  const int64_t rank_size = base_block + (rank < remainder ? 1 : 0);

  int thread_count = ppc::util::GetNumThreads();
  thread_count = std::max(1, std::min(thread_count, static_cast<int>(std::max<int64_t>(1, rank_size))));

  std::vector<std::thread> threads(thread_count);
  std::vector<double> thread_sums(thread_count, 0.0);

  const int64_t thread_block = rank_size / thread_count;
  const int64_t thread_remainder = rank_size % thread_count;

  int64_t current_begin = rank_begin;

  for (int thread_id = 0; thread_id < thread_count; ++thread_id) {
    const int64_t current_block = thread_block + (thread_id < thread_remainder ? 1 : 0);
    const int64_t current_end = current_begin + current_block;

    threads[thread_id] = std::thread([&, thread_id, current_begin, current_end]() {
      thread_sums[thread_id] = CalculateRange(current_begin, current_end, n, dimensions, h);
    });

    current_begin = current_end;
  }

  for (auto &thread : threads) {
    thread.join();
  }

  double local_sum = 0.0;

#pragma omp parallel for default(none) shared(thread_sums, thread_count) reduction(+ : local_sum) \
    num_threads(thread_count)
  for (int i = 0; i < thread_count; ++i) {
    local_sum += thread_sums[i];
  }

  double global_sum = 0.0;
  MPI_Allreduce(&local_sum, &global_sum, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

  GetOutput() = global_sum;

  for (int i = 0; i < dimensions; ++i) {
    GetOutput() *= h;
  }

  return true;
}

bool TelnovAIntegralRectangleALL::PostProcessingImpl() {
  return true;
}

}  // namespace telnov_a_integral_rectangle
