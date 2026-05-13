#include "kichanova_k_lin_system_by_conjug_grad/all/include/ops_all.hpp"

#include <omp.h>
#include <tbb/tbb.h>

#include <cmath>
#include <cstddef>
#include <vector>

#include "kichanova_k_lin_system_by_conjug_grad/common/include/common.hpp"

namespace kichanova_k_lin_system_by_conjug_grad {

namespace {

double DotProductHybrid(const std::vector<double> &a, const std::vector<double> &b, int n) {
  if (n <= 0) {
    return 0.0;
  }

  return tbb::parallel_reduce(tbb::blocked_range<int>(0, n, 1024), 0.0,
                              [&](const tbb::blocked_range<int> &range, double local_sum) {
#pragma omp simd reduction(+ : local_sum)
    for (int i = range.begin(); i < range.end(); ++i) {
      local_sum += a[i] * b[i];
    }
    return local_sum;
  }, [](double x, double y) { return x + y; });
}

void MatrixVectorProductHybrid(const std::vector<double> &a, const std::vector<double> &v, std::vector<double> &result,
                               int n) {
  if (n <= 0) {
    return;
  }

  const auto size_t stride = static_cast<size_t>(n);

  tbb::parallel_for(tbb::blocked_range<int>(0, n, 64), [&](const tbb::blocked_range<int> &range) {
    for (int i = range.begin(); i < range.end(); ++i) {
      const double *row = &a[i * stride];
      double sum = 0.0;

      int j = 0;
      for (; j <= n - 8; j += 8) {
        sum += row[j] * v[j];
        sum += row[j + 1] * v[j + 1];
        sum += row[j + 2] * v[j + 2];
        sum += row[j + 3] * v[j + 3];
        sum += row[j + 4] * v[j + 4];
        sum += row[j + 5] * v[j + 5];
        sum += row[j + 6] * v[j + 6];
        sum += row[j + 7] * v[j + 7];
      }
      for (; j < n; ++j) {
        sum += row[j] * v[j];
      }
      result[i] = sum;
    }
  });
}

void UpdateSolutionAndResidualHybrid(std::vector<double> &x, std::vector<double> &r, const std::vector<double> &p,
                                     const std::vector<double> &ap, double alpha, int n) {
  if (n <= 0) {
    return;
  }

  tbb::parallel_for(tbb::blocked_range<int>(0, n, 1024), [&](const tbb::blocked_range<int> &range) {
#pragma omp simd
    for (int i = range.begin(); i < range.end(); ++i) {
      x[i] += alpha * p[i];
      r[i] -= alpha * ap[i];
    }
  });
}

void UpdateDirectionHybrid(std::vector<double> &p, const std::vector<double> &r, double beta, int n) {
  if (n <= 0) {
    return;
  }

  tbb::parallel_for(tbb::blocked_range<int>(0, n, 1024), [&](const tbb::blocked_range<int> &range) {
#pragma omp simd
    for (int i = range.begin(); i < range.end(); ++i) {
      p[i] = r[i] + (beta * p[i]);
    }
  });
}

void InitializeVectorsHybrid(std::vector<double> &r, std::vector<double> &p, const std::vector<double> &b, int n) {
  tbb::parallel_for(tbb::blocked_range<int>(0, n, 1024), [&](const tbb::blocked_range<int> &range) {
    for (int i = range.begin(); i < range.end(); ++i) {
      r[i] = b[i];
      p[i] = r[i];
    }
  });
}

}  // namespace

KichanovaKLinSystemByConjugGradALL::KichanovaKLinSystemByConjugGradALL(const InType &in) {
  SetTypeOfTask(GetStaticTypeOfTask());
  GetInput() = in;
  GetOutput() = OutType();
}

bool KichanovaKLinSystemByConjugGradALL::ValidationImpl() {
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

bool KichanovaKLinSystemByConjugGradALL::PreProcessingImpl() {
  GetOutput().assign(GetInput().n, 0.0);
  return true;
}

bool KichanovaKLinSystemByConjugGradALL::RunImpl() {
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

  InitializeVectorsHybrid(r, p, b, n);

  double rr_old = DotProductHybrid(r, r, n);
  double residual_norm = std::sqrt(rr_old);

  if (residual_norm < epsilon) {
    return true;
  }

  int max_iter = n * 1000;

  for (int iter = 0; iter < max_iter; ++iter) {
    MatrixVectorProductHybrid(a, p, ap, n);

    double p_ap = DotProductHybrid(p, ap, n);
    if (std::abs(p_ap) < 1e-30) {
      break;
    }

    double alpha = rr_old / p_ap;
    UpdateSolutionAndResidualHybrid(x, r, p, ap, alpha, n);

    double rr_new = DotProductHybrid(r, r, n);
    residual_norm = std::sqrt(rr_new);

    if (residual_norm < epsilon) {
      break;
    }

    double beta = rr_new / rr_old;
    UpdateDirectionHybrid(p, r, beta, n);

    rr_old = rr_new;
  }

  return true;
}

bool KichanovaKLinSystemByConjugGradALL::PostProcessingImpl() {
  return true;
}

}  // namespace kichanova_k_lin_system_by_conjug_grad
