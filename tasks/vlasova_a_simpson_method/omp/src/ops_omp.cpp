#include "vlasova_a_simpson_method/omp/include/ops_omp.hpp"

#include <cmath>
#include <vector>
#include <omp.h>

#include "vlasova_a_simpson_method/common/include/common.hpp"

namespace vlasova_a_simpson_method {

VlasovaASimpsonMethodOMP::VlasovaASimpsonMethodOMP(InType in) : task_data_(std::move(in)) {
  SetTypeOfTask(GetStaticTypeOfTask());
  GetOutput() = 0.0;
}

bool VlasovaASimpsonMethodOMP::ValidationImpl() {
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

bool VlasovaASimpsonMethodOMP::PreProcessingImpl() {
  result_ = 0.0;
  GetOutput() = 0.0;

  size_t dim = task_data_.a.size();
  h_.resize(dim);
  dimensions_.resize(dim);
  total_points_ = 1;

  for (size_t i = 0; i < dim; ++i) {
    h_[i] = (task_data_.b[i] - task_data_.a[i]) / task_data_.n[i];
    dimensions_[i] = task_data_.n[i] + 1;
    total_points_ *= dimensions_[i];
  }

  weights_.resize(dim);
  for (size_t i = 0; i < dim; ++i) {
    int steps = task_data_.n[i];
    weights_[i].resize(steps + 1);
    for (int j = 0; j <= steps; ++j) {
      if (j == 0 || j == steps) {
        weights_[i][j] = 1.0;
      } else if (j % 2 == 0) {
        weights_[i][j] = 2.0;
      } else {
        weights_[i][j] = 4.0;
      }
    }
  }

  return true;
}

bool VlasovaASimpsonMethodOMP::RunImpl() {
  size_t dim = task_data_.a.size();
  double sum = 0.0;
  
  if (dim == 1) {
    int n0 = task_data_.n[0];
    double a0 = task_data_.a[0];
    double h0 = h_[0];
    const auto& w0 = weights_[0];
    const auto& func = task_data_.func;
    
    #pragma omp parallel for reduction(+:sum) schedule(static)
    for (int i = 0; i <= n0; ++i) {
      double x = a0 + i * h0;
      sum += w0[i] * func({x});
    }
    
    result_ = sum * h0 / 3.0;
    GetOutput() = result_;
    return true;
  }
  
  if (dim == 2) {
    int n0 = task_data_.n[0];
    int n1 = task_data_.n[1];
    double a0 = task_data_.a[0];
    double a1 = task_data_.a[1];
    double h0 = h_[0];
    double h1 = h_[1];
    const auto& w0 = weights_[0];
    const auto& w1 = weights_[1];
    const auto& func = task_data_.func;
    
    #pragma omp parallel for collapse(2) reduction(+:sum) schedule(static)
    for (int i = 0; i <= n0; ++i) {
      for (int j = 0; j <= n1; ++j) {
        double x = a0 + i * h0;
        double y = a1 + j * h1;
        sum += w0[i] * w1[j] * func({x, y});
      }
    }
    
    result_ = sum * h0 * h1 / 9.0;
    GetOutput() = result_;
    return true;
  }
  
  if (dim == 3) {
    int n0 = task_data_.n[0];
    int n1 = task_data_.n[1];
    int n2 = task_data_.n[2];
    double a0 = task_data_.a[0];
    double a1 = task_data_.a[1];
    double a2 = task_data_.a[2];
    double h0 = h_[0];
    double h1 = h_[1];
    double h2 = h_[2];
    const auto& w0 = weights_[0];
    const auto& w1 = weights_[1];
    const auto& w2 = weights_[2];
    const auto& func = task_data_.func;
    
    #pragma omp parallel for collapse(3) reduction(+:sum) schedule(static)
    for (int i = 0; i <= n0; ++i) {
      for (int j = 0; j <= n1; ++j) {
        for (int k = 0; k <= n2; ++k) {
          double x = a0 + i * h0;
          double y = a1 + j * h1;
          double z = a2 + k * h2;
          sum += w0[i] * w1[j] * w2[k] * func({x, y, z});
        }
      }
    }
    
    result_ = sum * h0 * h1 * h2 / 27.0;
    GetOutput() = result_;
    return true;
  }
  
  #pragma omp parallel reduction(+:sum)
  {
    std::vector<int> local_index(dim, 0);
    std::vector<double> local_point(dim);
    
    #pragma omp for schedule(static)
    for (size_t idx = 0; idx < total_points_; ++idx) {
      size_t temp_idx = idx;
      double weight = 1.0;
      
      for (size_t i = 0; i < dim; ++i) {
        int index_i = static_cast<int>(temp_idx % dimensions_[i]);
        temp_idx /= dimensions_[i];
        local_index[i] = index_i;
        weight *= weights_[i][index_i];
      }
      
      for (size_t i = 0; i < dim; ++i) {
        local_point[i] = task_data_.a[i] + (local_index[i] * h_[i]);
      }
      
      sum += weight * task_data_.func(local_point);
    }
  }
  
  double factor = 1.0;
  for (size_t i = 0; i < dim; ++i) {
    factor *= h_[i] / 3.0;
  }
  
  result_ = sum * factor;
  GetOutput() = result_;
  
  return true;
}

bool VlasovaASimpsonMethodOMP::PostProcessingImpl() {
  return std::isfinite(GetOutput());
}

}  // namespace vlasova_a_simpson_method