#include "kichanova_k_lin_system_by_conjug_grad/seq/include/ops_seq.hpp"

#include <cmath>
#include <cstddef>
#include <vector>

#include "kichanova_k_lin_system_by_conjug_grad/common/include/common.hpp"

namespace kichanova_k_lin_system_by_conjug_grad {

namespace {
double ComputeDotProduct(const std::vector<double> &a, const std::vector<double> &b, int n) {
  double result = 0.0;
  for (int i = 0; i < n; ++i) {
    result += a[i] * b[i];
  }
  return result;
}

void ComputeMatrixVectorProduct(const std::vector<double> &a, const std::vector<double> &v, std::vector<double> &result,
                                int n) {
  const auto stride = static_cast<size_t>(n);
  for (int i = 0; i < n; ++i) {
    double sum = 0.0;
    const double *a_row = &a[static_cast<size_t>(i) * stride];
    for (int j = 0; j < n; ++j) {
      sum += a_row[j] * v[j];
    }
    result[i] = sum;
  }
}

void UpdateSolution(std::vector<double> &x, const std::vector<double> &p, double alpha, int n) {
  for (int i = 0; i < n; ++i) {
    x[i] += alpha * p[i];
  }
}

void UpdateResidual(std::vector<double> &r, const std::vector<double> &ap, double alpha, int n) {
  for (int i = 0; i < n; ++i) {
    r[i] -= alpha * ap[i];
  }
}

void UpdateSearchDirection(std::vector<double> &p, const std::vector<double> &r, double beta, int n) {
  for (int i = 0; i < n; ++i) {
    p[i] = r[i] + (beta * p[i]);
  }
}
}  // namespace

KichanovaKLinSystemByConjugGradSEQ::KichanovaKLinSystemByConjugGradSEQ(const InType &in) {
  SetTypeOfTask(GetStaticTypeOfTask());
  GetInput() = in;
  GetOutput() = OutType();
}

bool KichanovaKLinSystemByConjugGradSEQ::ValidationImpl() {
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

bool KichanovaKLinSystemByConjugGradSEQ::PreProcessingImpl() {
  GetOutput().assign(GetInput().n, 0.0);
  return true;
}

bool KichanovaKLinSystemByConjugGradSEQ::RunImpl() {
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

  for (int i = 0; i < n; i++) {
    r[i] = b[i];
    p[i] = r[i];
  }

  double rr_old = ComputeDotProduct(r, r, n);
  double residual_norm = std::sqrt(rr_old);
  if (residual_norm < epsilon) {
    return true;
  }

  int max_iter = n * 1000;
  for (int iter = 0; iter < max_iter; iter++) {
    ComputeMatrixVectorProduct(a, p, ap, n);

    double p_ap = ComputeDotProduct(p, ap, n);
    if (std::abs(p_ap) < 1e-30) {
      break;
    }

    double alpha = rr_old / p_ap;
    UpdateSolution(x, p, alpha, n);
    UpdateResidual(r, ap, alpha, n);

    double rr_new = ComputeDotProduct(r, r, n);
    residual_norm = std::sqrt(rr_new);
    if (residual_norm < epsilon) {
      break;
    }

    double beta = rr_new / rr_old;
    UpdateSearchDirection(p, r, beta, n);

    rr_old = rr_new;
  }

  return true;
}

bool KichanovaKLinSystemByConjugGradSEQ::PostProcessingImpl() {
  return true;
}

}  // namespace kichanova_k_lin_system_by_conjug_grad
