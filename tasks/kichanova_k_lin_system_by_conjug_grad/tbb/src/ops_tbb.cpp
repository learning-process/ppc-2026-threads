#include "kichanova_k_lin_system_by_conjug_grad/tbb/include/ops_tbb.hpp"

#include <tbb/tbb.h>

#include <cmath>
#include <cstddef>
#include <vector>

#include "kichanova_k_lin_system_by_conjug_grad/common/include/common.hpp"

namespace kichanova_k_lin_system_by_conjug_grad {

namespace {

void ComputeMatrixVectorAndUpdateTBB(const std::vector<double> &a, const std::vector<double> &v,
                                     std::vector<double> &result, int n, double alpha, std::vector<double> &x,
                                     const std::vector<double> &p) {
  const auto stride = static_cast<size_t>(n);

  tbb::parallel_for(tbb::blocked_range<int>(0, n, 64), [&](const tbb::blocked_range<int> &range) {
    for (int i = range.begin(); i < range.end(); ++i) {
      double sum = 0.0;
      const double *a_row = &a[static_cast<size_t>(i) * stride];
      for (int j = 0; j < n; ++j) {
        sum += a_row[j] * v[j];
      }
      result[i] = sum;
      x[i] += alpha * p[i];
    }
  });
}

double ComputeResidualNorm(const std::vector<double> &r, int n) {
  double rr = tbb::parallel_reduce(tbb::blocked_range<int>(0, n, 1024), 0.0,
                                   [&](const tbb::blocked_range<int> &range, double local_sum) {
    for (int i = range.begin(); i < range.end(); ++i) {
      local_sum += r[i] * r[i];
    }
    return local_sum;
  }, [](double x, double y) { return x + y; });
  return std::sqrt(rr);
}

double ComputeDotProduct(const std::vector<double> &a, const std::vector<double> &b, int n) {
  return tbb::parallel_reduce(tbb::blocked_range<int>(0, n, 1024), 0.0,
                              [&](const tbb::blocked_range<int> &range, double local_sum) {
    for (int i = range.begin(); i < range.end(); ++i) {
      local_sum += a[i] * b[i];
    }
    return local_sum;
  }, [](double x, double y) { return x + y; });
}

void InitializeVectors(std::vector<double> &r, std::vector<double> &p, const std::vector<double> &b, int n) {
  tbb::parallel_for(tbb::blocked_range<int>(0, n, 1024), [&](const tbb::blocked_range<int> &range) {
    for (int i = range.begin(); i < range.end(); ++i) {
      r[i] = b[i];
      p[i] = r[i];
    }
  });
}

void UpdateSolutionAndResidual(std::vector<double> &x, std::vector<double> &r, const std::vector<double> &p,
                               const std::vector<double> &ap, double alpha, int n) {
  tbb::parallel_for(tbb::blocked_range<int>(0, n, 1024), [&](const tbb::blocked_range<int> &range) {
    for (int i = range.begin(); i < range.end(); ++i) {
      x[i] += alpha * p[i];
      r[i] -= alpha * ap[i];
    }
  });
}

void UpdateDirection(std::vector<double> &p, const std::vector<double> &r, double beta, int n) {
  tbb::parallel_for(tbb::blocked_range<int>(0, n, 1024), [&](const tbb::blocked_range<int> &range) {
    for (int i = range.begin(); i < range.end(); ++i) {
      p[i] = r[i] + (beta * p[i]);
    }
  });
}

}  // namespace

KichanovaKLinSystemByConjugGradTBB::KichanovaKLinSystemByConjugGradTBB(const InType &in) {
  SetTypeOfTask(GetStaticTypeOfTask());
  GetInput() = in;
  GetOutput() = OutType();
}

bool KichanovaKLinSystemByConjugGradTBB::ValidationImpl() {
  const InType &input_data = GetInput();
  if (input_data.n <= 0) {
    return false;
  }
  if (input_data.A.size() != static_cast<size_t>(input_data.n) * static_cast<size_t>(input_data.n)) {
    return false;
  }
  if (input_data.b.size() != static_cast<size_t>(input_data.n)) {
    return false;
  }
  return true;
}

bool KichanovaKLinSystemByConjugGradTBB::PreProcessingImpl() {
  GetOutput().assign(GetInput().n, 0.0);
  return true;
}

bool KichanovaKLinSystemByConjugGradTBB::RunImpl() {
  const InType &input_data = GetInput();
  OutType &x = GetOutput();

  int n = input_data.n;
  if (n == 0) {
    return false;
  }

  const std::vector<double> &a = input_data.A;
  const std::vector<double> &b = input_data.b;
  double epsilon = input_data.epsilon;

  std::vector<double> r(n);
  std::vector<double> p(n);
  std::vector<double> ap(n);

  InitializeVectors(r, p, b, n);

  double residual_norm = ComputeResidualNorm(r, n);
  if (residual_norm < epsilon) {
    return true;
  }

  double rr_old = residual_norm * residual_norm;
  int max_iter = n * 1000;

  for (int iter = 0; iter < max_iter; ++iter) {
    ComputeMatrixVectorAndUpdateTBB(a, p, ap, n, 0.0, x, p);

    double p_ap = ComputeDotProduct(p, ap, n);
    if (std::abs(p_ap) < 1e-30) {
      break;
    }

    double alpha = rr_old / p_ap;
    UpdateSolutionAndResidual(x, r, p, ap, alpha, n);

    residual_norm = ComputeResidualNorm(r, n);
    if (residual_norm < epsilon) {
      break;
    }

    double rr_new = residual_norm * residual_norm;
    double beta = rr_new / rr_old;
    UpdateDirection(p, r, beta, n);

    rr_old = rr_new;
  }

  return true;
}

bool KichanovaKLinSystemByConjugGradTBB::PostProcessingImpl() {
  return true;
}

}  // namespace kichanova_k_lin_system_by_conjug_grad
