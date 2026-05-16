#include "redkina_a_integral_simpson/all/include/ops_all.hpp"

#include <mpi.h>

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <functional>
#include <future>
#include <thread>
#include <vector>

#include "redkina_a_integral_simpson/common/include/common.hpp"

namespace redkina_a_integral_simpson {
namespace {

std::vector<std::vector<double>> PrecomputeWeights(const std::vector<int> &n) {
  const size_t dim = n.size();
  std::vector<std::vector<double>> weights(dim);
  for (size_t i = 0; i < dim; ++i) {
    const int ni = n[i];
    weights[i].resize(ni + 1);
    for (int idx = 0; idx <= ni; ++idx) {
      if (idx == 0 || idx == ni) {
        weights[i][idx] = 1.0;
      } else if (idx % 2 == 1) {
        weights[i][idx] = 4.0;
      } else {
        weights[i][idx] = 2.0;
      }
    }
  }
  return weights;
}

std::vector<size_t> ComputeStrides(const std::vector<int> &n) {
  const size_t dim = n.size();
  std::vector<size_t> strides(dim);
  if (dim == 0) {
    return strides;
  }
  strides[dim - 1] = 1;
  for (size_t i = dim - 1; i > 0; --i) {
    strides[i - 1] = strides[i] * static_cast<size_t>(n[i] + 1);
  }
  return strides;
}

double ComputeRangeSum(size_t start, size_t end, const std::vector<double> &a, const std::vector<double> &h,
                       const std::vector<std::vector<double>> &weights, const std::vector<size_t> &strides,
                       const std::function<double(const std::vector<double> &)> &func, size_t dim) {
  double sum = 0.0;
  std::vector<int> indices(dim);
  std::vector<double> point(dim);
  for (size_t idx = start; idx < end; ++idx) {
    size_t remainder = idx;
    for (size_t dim_idx = 0; dim_idx < dim; ++dim_idx) {
      indices[dim_idx] = static_cast<int>(remainder / strides[dim_idx]);
      remainder %= strides[dim_idx];
    }
    double w_prod = 1.0;
    for (size_t dim_idx = 0; dim_idx < dim; ++dim_idx) {
      const int i_idx = indices[dim_idx];
      point[dim_idx] = a[dim_idx] + (static_cast<double>(i_idx) * h[dim_idx]);
      w_prod *= weights[dim_idx][i_idx];
    }
    sum += w_prod * func(point);
  }
  return sum;
}

double ComputeLocalSumMPI(size_t local_start, size_t local_end, const std::vector<double> &a,
                          const std::vector<double> &h, const std::vector<std::vector<double>> &weights,
                          const std::vector<size_t> &strides,
                          const std::function<double(const std::vector<double> &)> &func, size_t dim) {
  const size_t local_size = local_end - local_start;
  if (local_size == 0) {
    return 0.0;
  }

  unsigned int hardware_threads = std::thread::hardware_concurrency();
  if (hardware_threads == 0) {
    hardware_threads = 2;
  }
  unsigned int num_threads = std::min(hardware_threads, static_cast<unsigned int>(local_size));

  if (num_threads == 1) {
    return ComputeRangeSum(local_start, local_end, a, h, weights, strides, func, dim);
  }

  std::vector<std::future<double>> futures;
  const size_t block_size = local_size / num_threads;
  const size_t rem_blocks = local_size % num_threads;
  size_t current_start = local_start;

  for (unsigned int thread_idx = 0; thread_idx < num_threads; ++thread_idx) {
    const size_t block_end = current_start + block_size + (thread_idx < rem_blocks ? 1 : 0);
    futures.push_back(
        std::async(std::launch::async, [&a, &h, &weights, &strides, &func, dim, current_start, block_end]() {
      return ComputeRangeSum(current_start, block_end, a, h, weights, strides, func, dim);
    }));
    current_start = block_end;
  }

  double total = 0.0;
  for (auto &f : futures) {
    total += f.get();
  }
  return total;
}

}  // namespace

RedkinaAIntegralSimpsonALL::RedkinaAIntegralSimpsonALL(const InType &in) {
  SetTypeOfTask(GetStaticTypeOfTask());
  GetInput() = in;
}

bool RedkinaAIntegralSimpsonALL::ValidationImpl() {
  const auto &in = GetInput();
  const size_t dim = in.a.size();

  if (dim == 0 || in.b.size() != dim || in.n.size() != dim) {
    return false;
  }
  for (size_t i = 0; i < dim; ++i) {
    if (in.a[i] >= in.b[i]) {
      return false;
    }
    if (in.n[i] <= 0 || in.n[i] % 2 != 0) {
      return false;
    }
  }
  return static_cast<bool>(in.func);
}

bool RedkinaAIntegralSimpsonALL::PreProcessingImpl() {
  const auto &in = GetInput();
  func_ = in.func;
  a_ = in.a;
  b_ = in.b;
  n_ = in.n;
  result_ = 0.0;
  return true;
}

bool RedkinaAIntegralSimpsonALL::RunImpl() {
  if (!func_) {
    return false;
  }
  const size_t dim = a_.size();
  if (dim == 0) {
    return false;
  }

  std::vector<double> h(dim);
  double h_prod = 1.0;
  for (size_t i = 0; i < dim; ++i) {
    h[i] = (b_[i] - a_[i]) / static_cast<double>(n_[i]);
    h_prod *= h[i];
  }

  const auto weights = PrecomputeWeights(n_);
  const auto strides = ComputeStrides(n_);
  if (strides.empty()) {
    return false;
  }

  const size_t total_points = strides[0] * static_cast<size_t>(n_[0] + 1);

  int rank = 0;
  int world_size = 1;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &world_size);

  const auto rank_u = static_cast<size_t>(rank);
  const auto size_u = static_cast<size_t>(world_size);
  const size_t base = total_points / size_u;
  const size_t rem = total_points % size_u;
  const size_t local_start = (rank_u * base) + std::min(rank_u, rem);
  const size_t local_end = local_start + base + (rank_u < rem ? 1 : 0);

  const double local_sum = ComputeLocalSumMPI(local_start, local_end, a_, h, weights, strides, func_, dim);

  double global_sum = 0.0;
  MPI_Allreduce(&local_sum, &global_sum, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

  double denominator = 1.0;
  for (size_t i = 0; i < dim; ++i) {
    denominator *= 3.0;
  }
  result_ = (h_prod / denominator) * global_sum;
  return true;
}

bool RedkinaAIntegralSimpsonALL::PostProcessingImpl() {
  GetOutput() = result_;
  return true;
}

}  // namespace redkina_a_integral_simpson
