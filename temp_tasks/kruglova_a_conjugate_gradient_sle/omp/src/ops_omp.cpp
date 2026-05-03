#include "kruglova_a_conjugate_gradient_sle/omp/include/ops_omp.hpp"

#include <cmath>
#include <cstddef>
#include <vector>

#include "kruglova_a_conjugate_gradient_sle/common/include/common.hpp"

namespace kruglova_a_conjugate_gradient_sle {

namespace {
void MatrixVectorMultiply(const std::vector<double> &a, const std::vector<double> &p, std::vector<double> &ap, int n) {
#pragma omp parallel for default(none) shared(a, p, ap, n)
  for (int i = 0; i < n; ++i) {
    ap[i] = 0.0;
    for (int j = 0; j < n; ++j) {
      ap[i] += a[(static_cast<size_t>(i) * n) + j] * p[j];
    }
  }
}
}  // namespace

KruglovaAConjGradSleOMP::KruglovaAConjGradSleOMP(const InType &in) {
  SetTypeOfTask(GetStaticTypeOfTask());
  GetInput() = in;
}

bool KruglovaAConjGradSleOMP::ValidationImpl() {
  const auto &in = GetInput();
  if (in.size <= 0) {
    return false;
  }
  if (in.A.size() != static_cast<size_t>(in.size) * static_cast<size_t>(in.size)) {
    return false;
  }
  if (in.b.size() != static_cast<size_t>(in.size)) {
    return false;
  }
  return true;
}

bool KruglovaAConjGradSleOMP::PreProcessingImpl() {
  GetOutput().assign(GetInput().size, 0.0);
  return true;
}

bool KruglovaAConjGradSleOMP::RunImpl() {
  const auto &a = GetInput().A;
  const auto &b = GetInput().b;
  int n = GetInput().size;
  auto &x = GetOutput();

  std::vector<double> r = b;
  std::vector<double> p = r;
  std::vector<double> ap(n, 0.0);

  double rsold = 0.0;
#pragma omp parallel for default(none) shared(r, n) reduction(+ : rsold)
  for (int i = 0; i < n; ++i) {
    rsold += r[i] * r[i];
  }

  const double tolerance = 1e-8;

  for (int iter = 0; iter < n * 2; ++iter) {
    MatrixVectorMultiply(a, p, ap, n);

    double p_ap = 0.0;
#pragma omp parallel for default(none) shared(p, ap, n) reduction(+ : p_ap)
    for (int i = 0; i < n; ++i) {
      p_ap += p[i] * ap[i];
    }

    if (std::abs(p_ap) < 1e-15) {
      break;
    }

    double alpha = rsold / p_ap;

#pragma omp parallel for default(none) shared(x, r, ap, p, alpha, n)
    for (int i = 0; i < n; ++i) {
      x[i] += alpha * p[i];
      r[i] -= alpha * ap[i];
    }

    double rsnew = 0.0;
#pragma omp parallel for default(none) shared(r, n) reduction(+ : rsnew)
    for (int i = 0; i < n; ++i) {
      rsnew += r[i] * r[i];
    }

    if (std::sqrt(rsnew) < tolerance) {
      break;
    }

    double beta = rsnew / rsold;
#pragma omp parallel for default(none) shared(p, r, beta, n)
    for (int i = 0; i < n; ++i) {
      p[i] = r[i] + (beta * p[i]);
    }
    rsold = rsnew;
  }
  return true;
}

bool KruglovaAConjGradSleOMP::PostProcessingImpl() {
  return true;
}

}  // namespace kruglova_a_conjugate_gradient_sle
