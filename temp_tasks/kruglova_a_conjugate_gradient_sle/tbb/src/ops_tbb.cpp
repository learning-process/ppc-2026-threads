#include "kruglova_a_conjugate_gradient_sle/tbb/include/ops_tbb.hpp"

#include <oneapi/tbb/blocked_range.h>
#include <oneapi/tbb/parallel_for.h>
#include <oneapi/tbb/parallel_reduce.h>

#include <cmath>
#include <cstddef>
#include <vector>

#include "kruglova_a_conjugate_gradient_sle/common/include/common.hpp"

namespace kruglova_a_conjugate_gradient_sle {

namespace {
void MatrixVectorMultiply(const std::vector<double> &a, const std::vector<double> &p, std::vector<double> &ap, int n) {
  tbb::parallel_for(tbb::blocked_range<int>(0, n, 256), [&](const tbb::blocked_range<int> &range) {
    for (int i = range.begin(); i < range.end(); ++i) {
      double sum = 0.0;
      const size_t row_offset = static_cast<size_t>(i) * n;
      for (int j = 0; j < n; ++j) {
        sum += a[row_offset + j] * p[j];
      }
      ap[i] = sum;
    }
  });
}

double DotProduct(const std::vector<double> &v1, const std::vector<double> &v2, int n) {
  return tbb::parallel_reduce(tbb::blocked_range<int>(0, n, 512), 0.0,
                              [&](const tbb::blocked_range<int> &range, double init) {
    for (int i = range.begin(); i < range.end(); ++i) {
      init += v1[i] * v2[i];
    }
    return init;
  }, [](double a, double b) { return a + b; });
}
}  // namespace

KruglovaAConjGradSleTBB::KruglovaAConjGradSleTBB(const InType &in) {
  SetTypeOfTask(GetStaticTypeOfTask());
  GetInput() = in;
}

bool KruglovaAConjGradSleTBB::ValidationImpl() {
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

bool KruglovaAConjGradSleTBB::PreProcessingImpl() {
  GetOutput().assign(GetInput().size, 0.0);
  return true;
}

bool KruglovaAConjGradSleTBB::RunImpl() {
  const auto &a = GetInput().A;
  const auto &b = GetInput().b;
  int n = GetInput().size;
  auto &x = GetOutput();

  std::vector<double> r = b;
  std::vector<double> p = r;
  std::vector<double> ap(n, 0.0);

  double rsold = DotProduct(r, r, n);

  const double tolerance = 1e-8;

  for (int iter = 0; iter < n * 2; ++iter) {
    MatrixVectorMultiply(a, p, ap, n);

    double p_ap = DotProduct(p, ap, n);

    if (std::abs(p_ap) < 1e-15) {
      break;
    }

    double alpha = rsold / p_ap;

    tbb::parallel_for(tbb::blocked_range<int>(0, n, 1024), [&](const tbb::blocked_range<int> &range) {
      for (int i = range.begin(); i < range.end(); ++i) {
        x[i] += alpha * p[i];
        r[i] -= alpha * ap[i];
      }
    });

    double rsnew = DotProduct(r, r, n);

    if (std::sqrt(rsnew) < tolerance) {
      break;
    }

    double beta = rsnew / rsold;

    tbb::parallel_for(tbb::blocked_range<int>(0, n, 1024), [&](const tbb::blocked_range<int> &range) {
      for (int i = range.begin(); i < range.end(); ++i) {
        p[i] = r[i] + (beta * p[i]);
      }
    });

    rsold = rsnew;
  }
  return true;
}

bool KruglovaAConjGradSleTBB::PostProcessingImpl() {
  return true;
}

}  // namespace kruglova_a_conjugate_gradient_sle
