#include "smyshlaev_a_sle_cg_seq/omp/include/ops_omp.hpp"

#include <cmath>
#include <cstddef>
#include <vector>

#include "smyshlaev_a_sle_cg_seq/common/include/common.hpp"
#include "util/include/util.hpp"

namespace smyshlaev_a_sle_cg_seq {

namespace {

double ComputeDotProduct(const std::vector<double> &v1, const std::vector<double> &v2, int num_threads) {
  double result = 0.0;
  int n = static_cast<int>(v1.size());
#pragma omp parallel for default(none) shared(n, v1, v2) num_threads(num_threads) reduction(+ : result)
  for (int i = 0; i < n; ++i) {
    result += v1[i] * v2[i];
  }
  return result;
}
void ComputeAp(const std::vector<double> &matrix, const std::vector<double> &p, std::vector<double> &ap, int n,
               int num_threads) {
#pragma omp parallel for default(none) shared(n, matrix, p, ap) num_threads(num_threads)
  for (int i = 0; i < n; ++i) {
    double sum = 0.0;
    for (int j = 0; j < n; ++j) {
      sum += matrix[(i * n) + j] * p[j];
    }
    ap[i] = sum;
  }
}
double UpdateResultAndResidual(std::vector<double> &result, std::vector<double> &r, const std::vector<double> &p,
                               const std::vector<double> &ap, double alpha, int num_threads) {
  double rs_new = 0.0;
  int n = static_cast<int>(result.size());
#pragma omp parallel for default(none) shared(n, result, r, p, ap, alpha) num_threads(num_threads) reduction(+ : rs_new)
  for (int i = 0; i < n; ++i) {
    result[i] += alpha * p[i];
    r[i] -= alpha * ap[i];
    rs_new += r[i] * r[i];
  }
  return rs_new;
}

void UpdateP(std::vector<double> &p, const std::vector<double> &r, double beta, int num_threads) {
  int n = static_cast<int>(p.size());
#pragma omp parallel for default(none) shared(n, p, r, beta) num_threads(num_threads)
  for (int i = 0; i < n; ++i) {
    p[i] = r[i] + (beta * p[i]);
  }
}

double ComputeDotProduct(const std::vector<double> &v1, const std::vector<double> &v2) {
  double result = 0.0;
  int n = static_cast<int>(v1.size());
  for (int i = 0; i < n; ++i) {
    result += v1[i] * v2[i];
  }
  return result;
}
void ComputeAp(const std::vector<double> &matrix, const std::vector<double> &p, std::vector<double> &ap, int n) {
  for (int i = 0; i < n; ++i) {
    double sum = 0.0;
    for (int j = 0; j < n; ++j) {
      sum += matrix[(i * n) + j] * p[j];
    }
    ap[i] = sum;
  }
}
double UpdateResultAndResidual(std::vector<double> &result, std::vector<double> &r, const std::vector<double> &p,
                               const std::vector<double> &ap, double alpha) {
  double rs_new = 0.0;
  int n = static_cast<int>(result.size());
  for (int i = 0; i < n; ++i) {
    result[i] += alpha * p[i];
    r[i] -= alpha * ap[i];
    rs_new += r[i] * r[i];
  }
  return rs_new;
}

void UpdateP(std::vector<double> &p, const std::vector<double> &r, double beta) {
  int n = static_cast<int>(p.size());
  for (int i = 0; i < n; ++i) {
    p[i] = r[i] + (beta * p[i]);
  }
}

}  // namespace

SmyshlaevASleCgTaskOMP::SmyshlaevASleCgTaskOMP(const InType &in) {
  SetTypeOfTask(GetStaticTypeOfTask());
  GetInput() = in;
}

bool SmyshlaevASleCgTaskOMP::ValidationImpl() {
  const auto &a = GetInput().A;
  const auto &b = GetInput().b;
  if (a.empty() || b.empty()) {
    return false;
  }
  if (a.size() != b.size()) {
    return false;
  }
  if (a.size() != a[0].size()) {
    return false;
  }
  return true;
}

bool SmyshlaevASleCgTaskOMP::PreProcessingImpl() {
  const auto &a = GetInput().A;
  size_t n = a.size();
  flat_A_.resize(n * n);
  for (size_t i = 0; i < n; ++i) {
    for (size_t j = 0; j < n; ++j) {
      flat_A_[(i * n) + j] = a[i][j];
    }
  }
  return true;
}

bool SmyshlaevASleCgTaskOMP::RunSequential() {
  const auto &b = GetInput().b;
  int n = static_cast<int>(b.size());
  std::vector<double> r = b;
  std::vector<double> p = r;
  std::vector<double> ap(n, 0.0);
  std::vector<double> result(n, 0.0);

  double rs_old = 0.0;
  for (int i = 0; i < n; ++i) {
    rs_old += r[i] * r[i];
  }

  const double epsilon = 1e-9;
  if (std::sqrt(rs_old) < epsilon) {
    GetOutput() = result;
    return true;
  }

  const int max_iterations = n * 2;
  for (int iter = 0; iter < max_iterations; ++iter) {
    ComputeAp(flat_A_, p, ap, n);
    double p_ap = ComputeDotProduct(p, ap);
    if (std::abs(p_ap) < 1e-15) {
      break;
    }

    double alpha = rs_old / p_ap;
    double rs_new = UpdateResultAndResidual(result, r, p, ap, alpha);
    if (std::sqrt(rs_new) < epsilon) {
      break;
    }

    double beta = rs_new / rs_old;
    UpdateP(p, r, beta);
    rs_old = rs_new;
  }

  GetOutput() = result;
  return true;
}

bool SmyshlaevASleCgTaskOMP::RunParallel(int num_threads) {
  const auto &b = GetInput().b;
  int n = static_cast<int>(b.size());
  std::vector<double> r = b;
  std::vector<double> p = r;
  std::vector<double> ap(n, 0.0);
  std::vector<double> result(n, 0.0);

  double rs_old = 0.0;

#pragma omp parallel for default(none) shared(n, r) num_threads(num_threads) reduction(+ : rs_old)
  for (int i = 0; i < n; ++i) {
    rs_old += r[i] * r[i];
  }

  const double epsilon = 1e-9;
  if (std::sqrt(rs_old) < epsilon) {
    GetOutput() = result;
    return true;
  }

  const int max_iterations = n * 2;
  for (int iter = 0; iter < max_iterations; ++iter) {
    ComputeAp(flat_A_, p, ap, n, num_threads);
    double p_ap = ComputeDotProduct(p, ap, num_threads);
    if (std::abs(p_ap) < 1e-15) {
      break;
    }

    double alpha = rs_old / p_ap;
    double rs_new = UpdateResultAndResidual(result, r, p, ap, alpha, num_threads);
    if (std::sqrt(rs_new) < epsilon) {
      break;
    }

    double beta = rs_new / rs_old;
    UpdateP(p, r, beta, num_threads);
    rs_old = rs_new;
  }

  GetOutput() = result;
  return true;
}

bool SmyshlaevASleCgTaskOMP::RunImpl() {
  const auto &b = GetInput().b;
  if (b.empty()) {
    return true;
  }
  int n = static_cast<int>(b.size());

  int num_threads = ppc::util::GetNumThreads();
  if (n < num_threads * 100) {
    num_threads = 1;
  }
  if (num_threads == 1) {
    return RunSequential();
  }

  return RunParallel(num_threads);
}

bool SmyshlaevASleCgTaskOMP::PostProcessingImpl() {
  return true;
}

}  // namespace smyshlaev_a_sle_cg_seq
