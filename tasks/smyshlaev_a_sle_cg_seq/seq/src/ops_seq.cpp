#include "smyshlaev_a_sle_cg_seq/seq/include/ops_seq.hpp"

#include <cmath>
#include <cstddef>
#include <vector>

#include "smyshlaev_a_sle_cg_seq/common/include/common.hpp"

namespace smyshlaev_a_sle_cg_seq {

namespace {

double ComputeDotProduct(const std::vector<double> &v1, const std::vector<double> &v2) {
  double result = 0.0;
  size_t n = v1.size();
  for (size_t i = 0; i < n; ++i) {
    result += v1[i] * v2[i];
  }
  return result;
}

void ComputeAp(const std::vector<double> &matrix, const std::vector<double> &p, std::vector<double> &ap, size_t n) {
  for (size_t i = 0; i < n; ++i) {
    double sum = 0.0;
    for (size_t j = 0; j < n; ++j) {
      sum += matrix[(i * n) + j] * p[j];
    }
    ap[i] = sum;
  }
}

double UpdateResultAndResidual(std::vector<double> &result, std::vector<double> &r, const std::vector<double> &p,
                               const std::vector<double> &ap, double alpha) {
  double rs_new = 0.0;
  size_t n = result.size();
  for (size_t i = 0; i < n; ++i) {
    result[i] += alpha * p[i];
    r[i] -= alpha * ap[i];
    rs_new += r[i] * r[i];
  }
  return rs_new;
}

void UpdateP(std::vector<double> &p, const std::vector<double> &r, double beta) {
  size_t n = p.size();
  for (size_t i = 0; i < n; ++i) {
    p[i] = r[i] + (beta * p[i]);
  }
}

}  // namespace

SmyshlaevASleCgTaskSEQ::SmyshlaevASleCgTaskSEQ(const InType &in) {
  SetTypeOfTask(GetStaticTypeOfTask());
  GetInput() = in;
}

bool SmyshlaevASleCgTaskSEQ::ValidationImpl() {
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

bool SmyshlaevASleCgTaskSEQ::PreProcessingImpl() {
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

bool SmyshlaevASleCgTaskSEQ::RunImpl() {
  const auto &b = GetInput().b;
  size_t n = b.size();

  if (n == 0) {
    return true;
  }

  std::vector<double> r = b;
  std::vector<double> p = r;
  std::vector<double> ap(n, 0.0);
  std::vector<double> result(n, 0.0);

  double rs_old = ComputeDotProduct(r, r);

  const int max_iterations = static_cast<int>(n) * 2;
  const double epsilon = 1e-9;

  if (std::sqrt(rs_old) < epsilon) {
    GetOutput() = result;
    return true;
  }

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

bool SmyshlaevASleCgTaskSEQ::PostProcessingImpl() {
  return true;
}

}  // namespace smyshlaev_a_sle_cg_seq
