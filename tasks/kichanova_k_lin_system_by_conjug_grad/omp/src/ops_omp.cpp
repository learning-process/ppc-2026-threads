#include "kichanova_k_lin_system_by_conjug_grad/omp/include/ops_omp.hpp"

#include <cmath>
#include <cstddef>
#include <vector>

#include "kichanova_k_lin_system_by_conjug_grad/common/include/common.hpp"

namespace kichanova_k_lin_system_by_conjug_grad {

namespace {
double ComputeDotProductOMP(const std::vector<double> &a, const std::vector<double> &b, int n) {
  double result = 0.0;
#pragma omp parallel for default(none) shared(a, b, n) reduction(+ : result) schedule(static)
  for (int i = 0; i < n; ++i) {
    result += a[i] * b[i];
  }
  return result;
}

void ComputeMatrixVectorProductOMP(const std::vector<double> &a, const std::vector<double> &v,
                                   std::vector<double> &result, int n) {
  const auto stride = static_cast<size_t>(n);
#pragma omp parallel for default(none) shared(a, v, result, n, stride) schedule(static)
  for (int i = 0; i < n; ++i) {
    double sum = 0.0;
    const double *a_row = &a[static_cast<size_t>(i) * stride];
    for (int j = 0; j < n; ++j) {
      sum += a_row[j] * v[j];
    }
    result[i] = sum;
  }
}

void UpdateSolutionOMP(std::vector<double> &x, const std::vector<double> &p, double alpha, int n) {
#pragma omp parallel for default(none) shared(x, p, alpha, n) schedule(static)
  for (int i = 0; i < n; ++i) {
    x[i] += alpha * p[i];
  }
}

void UpdateResidualOMP(std::vector<double> &r, const std::vector<double> &ap, double alpha, int n) {
#pragma omp parallel for default(none) shared(r, ap, alpha, n) schedule(static)
  for (int i = 0; i < n; ++i) {
    r[i] -= alpha * ap[i];
  }
}

void UpdateSearchDirectionOMP(std::vector<double> &p, const std::vector<double> &r, double beta, int n) {
#pragma omp parallel for default(none) shared(p, r, beta, n) schedule(static)
  for (int i = 0; i < n; ++i) {
    p[i] = r[i] + (beta * p[i]);
  }
}
}  // namespace

KichanovaKLinSystemByConjugGradOMP::KichanovaKLinSystemByConjugGradOMP(const InType &in) {
  SetTypeOfTask(GetStaticTypeOfTask());
  GetInput() = in;
  GetOutput() = OutType();
}

bool KichanovaKLinSystemByConjugGradOMP::ValidationImpl() {
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

bool KichanovaKLinSystemByConjugGradOMP::PreProcessingImpl() {
  GetOutput().assign(GetInput().n, 0.0);
  return true;
}

bool KichanovaKLinSystemByConjugGradOMP::RunImpl() {
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

#pragma omp parallel for default(none) shared(r, b, p, n) schedule(static)
  for (int i = 0; i < n; i++) {
    r[i] = b[i];
    p[i] = r[i];
  }

  double rr_old = ComputeDotProductOMP(r, r, n);
  double residual_norm = std::sqrt(rr_old);
  if (residual_norm < epsilon) {
    return true;
  }

  int max_iter = n * 1000;

  for (int iter = 0; iter < max_iter; iter++) {
    ComputeMatrixVectorProductOMP(a, p, ap, n);

    double p_ap = ComputeDotProductOMP(p, ap, n);
    if (std::abs(p_ap) < 1e-30) {
      break;
    }

    double alpha = rr_old / p_ap;

    UpdateSolutionOMP(x, p, alpha, n);

    UpdateResidualOMP(r, ap, alpha, n);

    double rr_new = ComputeDotProductOMP(r, r, n);
    residual_norm = std::sqrt(rr_new);
    if (residual_norm < epsilon) {
      break;
    }

    double beta = rr_new / rr_old;

    UpdateSearchDirectionOMP(p, r, beta, n);

    rr_old = rr_new;
  }

  return true;
}

bool KichanovaKLinSystemByConjugGradOMP::PostProcessingImpl() {
  return true;
}

}  // namespace kichanova_k_lin_system_by_conjug_grad
