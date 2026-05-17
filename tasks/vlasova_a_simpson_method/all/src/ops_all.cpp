#include "vlasova_a_simpson_method/all/include/ops_all.hpp"

#include <mpi.h>
#include <tbb/blocked_range.h>
#include <tbb/parallel_reduce.h>

#include <cmath>
#include <cstddef>
#include <utility>
#include <vector>

#include "vlasova_a_simpson_method/common/include/common.hpp"

namespace vlasova_a_simpson_method {

VlasovaASimpsonMethodALL::VlasovaASimpsonMethodALL(InType in) : task_data_(std::move(in)) {
  SetTypeOfTask(GetStaticTypeOfTask());
  GetOutput() = 0.0;
}

bool VlasovaASimpsonMethodALL::ValidationImpl() {
  size_t dim = task_data_.a.size();

  if (dim == 0 || dim != task_data_.b.size() || dim != task_data_.n.size()) {
    return false;
  }

  for (size_t i = 0; i < dim; ++i) {
    if (task_data_.a[i] >= task_data_.b[i]) {
      return false;
    }
    if (task_data_.n[i] <= 0 || task_data_.n[i] % 2 != 0) {
      return false;
    }
  }

  if (!task_data_.func) {
    return false;
  }

  return GetOutput() == 0.0;
}

bool VlasovaASimpsonMethodALL::PreProcessingImpl() {
  result_ = 0.0;
  GetOutput() = 0.0;

  size_t dim = task_data_.a.size();
  h_.resize(dim);
  dimensions_.resize(dim);

  for (size_t i = 0; i < dim; ++i) {
    h_[i] = (task_data_.b[i] - task_data_.a[i]) / task_data_.n[i];
    dimensions_[i] = task_data_.n[i] + 1;
  }

  return true;
}

void VlasovaASimpsonMethodALL::ComputeWeight(const std::vector<int> &index, double &weight) const {
  weight = 1.0;
  size_t dim = index.size();

  for (size_t i = 0; i < dim; ++i) {
    int idx = index[i];
    int steps = task_data_.n[i];

    if (idx == 0 || idx == steps) {
      weight *= 1.0;
    } else if (idx % 2 == 0) {
      weight *= 2.0;
    } else {
      weight *= 4.0;
    }
  }
}

void VlasovaASimpsonMethodALL::ComputePoint(const std::vector<int> &index, std::vector<double> &point) const {
  size_t dim = index.size();
  point.resize(dim);

  for (size_t i = 0; i < dim; ++i) {
    point[i] = task_data_.a[i] + (index[i] * h_[i]);
  }
}

bool VlasovaASimpsonMethodALL::RunImpl() {
  int rank = 0;
  int size = 1;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &size);

  size_t dim = task_data_.a.size();

  size_t total_points = 1;
  for (size_t i = 0; i < dim; ++i) {
    total_points *= static_cast<size_t>(dimensions_[i]);
  }

  const auto urank = static_cast<size_t>(rank);
  const auto usize = static_cast<size_t>(size);
  const size_t chunk = total_points / usize;
  const size_t rem = total_points % usize;

  size_t begin = 0;
  size_t end = 0;
  if (urank < rem) {
    begin = urank * (chunk + 1);
    end = begin + (chunk + 1);
  } else {
    begin = (rem * (chunk + 1)) + ((urank - rem) * chunk);
    end = begin + chunk;
  }

  double local_sum = tbb::parallel_reduce(tbb::blocked_range<size_t>(begin, end), 0.0,
                                          [this, dim](const tbb::blocked_range<size_t> &range, double acc) {
    std::vector<int> cur_index(dim, 0);
    std::vector<double> cur_point;
    double local_weight = 0.0;

    for (size_t idx = range.begin(); idx != range.end(); ++idx) {
      size_t temp = idx;
      for (size_t i = 0; i < dim; ++i) {
        cur_index[i] = static_cast<int>(temp % static_cast<size_t>(dimensions_[i]));
        temp /= static_cast<size_t>(dimensions_[i]);
      }
      ComputeWeight(cur_index, local_weight);
      ComputePoint(cur_index, cur_point);
      acc += local_weight * task_data_.func(cur_point);
    }
    return acc;
  }, [](double x, double y) { return x + y; });

  double global_sum = 0.0;
  MPI_Allreduce(&local_sum, &global_sum, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

  double factor = 1.0;
  for (size_t i = 0; i < dim; ++i) {
    factor *= h_[i] / 3.0;
  }
  result_ = global_sum * factor;
  GetOutput() = result_;

  return true;
}

bool VlasovaASimpsonMethodALL::PostProcessingImpl() {
  return std::isfinite(GetOutput());
}

}  // namespace vlasova_a_simpson_method
