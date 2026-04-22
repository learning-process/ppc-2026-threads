#include "smyshlaev_a_sle_cg_seq/tbb/include/ops_tbb.hpp"  // Изменено на tbb

#include <tbb/tbb.h>

#include <cmath>
#include <cstddef>
#include <functional>
#include <vector>

#include "oneapi/tbb/parallel_for.h"
#include "smyshlaev_a_sle_cg_seq/common/include/common.hpp"
#include "util/include/util.hpp"

namespace smyshlaev_a_sle_cg_seq {

namespace {

double ComputeDotProductTbb(const std::vector<double> &v1, const std::vector<double> &v2) {
  int n = static_cast<int>(v1.size());
  return tbb::parallel_reduce(tbb::blocked_range<int>(0, n), 0.0,
                              [&](const tbb::blocked_range<int> &range, double init) -> double {
    for (int i = range.begin(); i != range.end(); ++i) {
      init += v1[i] * v2[i];
    }
    return init;
  }, std::plus<>());
}

void ComputeApTbb(const std::vector<double> &matrix, const std::vector<double> &p, std::vector<double> &ap, int n) {
  tbb::parallel_for(tbb::blocked_range<int>(0, n), [&](const tbb::blocked_range<int> &range) {
    for (int i = range.begin(); i != range.end(); ++i) {
      double sum = 0.0;
      for (int j = 0; j < n; ++j) {
        sum += matrix[(i * n) + j] * p[j];
      }
      ap[i] = sum;
    }
  });
}

double UpdateResultAndResidualTbb(std::vector<double> &result, std::vector<double> &r, const std::vector<double> &p,
                                  const std::vector<double> &ap, double alpha) {
  int n = static_cast<int>(result.size());
  return tbb::parallel_reduce(tbb::blocked_range<int>(0, n), 0.0,
                              [&](const tbb::blocked_range<int> &range, double init) -> double {
    for (int i = range.begin(); i != range.end(); ++i) {
      result[i] += alpha * p[i];
      r[i] -= alpha * ap[i];
      init += r[i] * r[i];
    }
    return init;
  }, std::plus<>());
}

void UpdatePTbb(std::vector<double> &p, const std::vector<double> &r, double beta) {
  int n = static_cast<int>(p.size());
  tbb::parallel_for(tbb::blocked_range<int>(0, n), [&](const tbb::blocked_range<int> &range) {
    for (int i = range.begin(); i != range.end(); ++i) {
      p[i] = r[i] + (beta * p[i]);
    }
  });
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

SmyshlaevASleCgTaskTBB::SmyshlaevASleCgTaskTBB(const InType &in) {
  SetTypeOfTask(GetStaticTypeOfTask());
  GetInput() = in;
}

bool SmyshlaevASleCgTaskTBB::ValidationImpl() {
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

bool SmyshlaevASleCgTaskTBB::PreProcessingImpl() {
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

bool SmyshlaevASleCgTaskTBB::RunSequential() {
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

bool SmyshlaevASleCgTaskTBB::RunParallel(int num_threads) {
  tbb::global_control thread_limit(tbb::global_control::max_allowed_parallelism, num_threads);

  const auto &b = GetInput().b;
  int n = static_cast<int>(b.size());
  std::vector<double> r = b;
  std::vector<double> p = r;
  std::vector<double> ap(n, 0.0);
  std::vector<double> result(n, 0.0);

  double rs_old = tbb::parallel_reduce(tbb::blocked_range<int>(0, n), 0.0,
                                       [&](const tbb::blocked_range<int> &range, double init) -> double {
    for (int i = range.begin(); i != range.end(); ++i) {
      init += r[i] * r[i];
    }
    return init;
  }, std::plus<>());

  const double epsilon = 1e-9;
  if (std::sqrt(rs_old) < epsilon) {
    GetOutput() = result;
    return true;
  }

  const int max_iterations = n * 2;
  for (int iter = 0; iter < max_iterations; ++iter) {
    ComputeApTbb(flat_A_, p, ap, n);
    double p_ap = ComputeDotProductTbb(p, ap);
    if (std::abs(p_ap) < 1e-15) {
      break;
    }

    double alpha = rs_old / p_ap;
    double rs_new = UpdateResultAndResidualTbb(result, r, p, ap, alpha);
    if (std::sqrt(rs_new) < epsilon) {
      break;
    }

    double beta = rs_new / rs_old;
    UpdatePTbb(p, r, beta);
    rs_old = rs_new;
  }

  GetOutput() = result;
  return true;
}

bool SmyshlaevASleCgTaskTBB::RunImpl() {
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

bool SmyshlaevASleCgTaskTBB::PostProcessingImpl() {
  return true;
}

}  // namespace smyshlaev_a_sle_cg_seq
