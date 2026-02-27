#include "kruglova_a_conjugate_gradient_sle/seq/include/ops_seq.hpp"

#include <cmath>
#include <vector>

namespace kruglova_a_conjugate_gradient_sle {

KruglovaAConjGradSleSEQ::KruglovaAConjGradSleSEQ(const InType &in) {
  SetTypeOfTask(GetStaticTypeOfTask());
  GetInput() = in;
}

bool KruglovaAConjGradSleSEQ::ValidationImpl() {
  const auto &in = GetInput();
  if (in.size <= 0) {
    return false;
  }

  if (in.A.size() != static_cast<size_t>(in.size * in.size)) {
    return false;
  }
  if (in.b.size() != static_cast<size_t>(in.size)) {
    return false;
  }
  return true;
}

bool KruglovaAConjGradSleSEQ::PreProcessingImpl() {
  GetOutput().assign(GetInput().size, 0.0);
  return true;
}

bool KruglovaAConjGradSleSEQ::RunImpl() {
  const auto &A = GetInput().A;
  const auto &b = GetInput().b;
  int n = GetInput().size;
  auto &x = GetOutput();

  std::vector<double> r = b;
  std::vector<double> p = r;
  std::vector<double> Ap(n, 0.0);

  double rsold = 0.0;
  for (int i = 0; i < n; ++i) {
    rsold += r[i] * r[i];
  }

  const double tolerance = 1e-8;

  for (int iter = 0; iter < n * 2; ++iter) {
    for (int i = 0; i < n; ++i) {
      Ap[i] = 0.0;
      for (int j = 0; j < n; ++j) {
        Ap[i] += A[i * n + j] * p[j];
      }
    }

    double pAp = 0.0;
    for (int i = 0; i < n; ++i) {
      pAp += p[i] * Ap[i];
    }

    if (pAp == 0.0) {
      break;
    }

    double alpha = rsold / pAp;

    for (int i = 0; i < n; ++i) {
      x[i] += alpha * p[i];
      r[i] -= alpha * Ap[i];
    }

    double rsnew = 0.0;
    for (int i = 0; i < n; ++i) {
      rsnew += r[i] * r[i];
    }

    if (std::sqrt(rsnew) < tolerance) {
      break;
    }

    for (int i = 0; i < n; ++i) {
      p[i] = r[i] + (rsnew / rsold) * p[i];
    }
    rsold = rsnew;
  }

  return true;
}

bool KruglovaAConjGradSleSEQ::PostProcessingImpl() {
  return true;
}

}  // namespace kruglova_a_conjugate_gradient_sle
