#include "kruglova_a_conjugate_gradient_sle/stl/include/ops_stl.hpp"

#include <algorithm>
#include <cmath>
#include <execution>
#include <numeric>
#include <vector>

#include "kruglova_a_conjugate_gradient_sle/common/include/common.hpp"

namespace kruglova_a_conjugate_gradient_sle {

// Параллельное скалярное произведение через Parallel STL
static double ParallelDotProduct(const std::vector<double> &v1, const std::vector<double> &v2) {
  if (v1.empty()) {
    return 0.0;
  }
  return std::transform_reduce(std::execution::par, v1.begin(), v1.end(), v2.begin(), 0.0);
}

KruglovaAConjGradSleSTL::KruglovaAConjGradSleSTL(const InType &in) {
  SetTypeOfTask(GetStaticTypeOfTask());
  GetInput() = in;
}

bool KruglovaAConjGradSleSTL::ValidationImpl() {
  const auto &in = GetInput();
  return in.size > 0 && in.A.size() == static_cast<size_t>(in.size) * in.size &&
         in.b.size() == static_cast<size_t>(in.size);
}

bool KruglovaAConjGradSleSTL::PreProcessingImpl() {
  GetOutput().assign(GetInput().size, 0.0);
  return true;
}

bool KruglovaAConjGradSleSTL::RunImpl() {
  const auto &a = GetInput().A;
  const auto &b = GetInput().b;
  const int n = GetInput().size;
  auto &x = GetOutput();

  std::vector<double> r = b;
  std::vector<double> p = r;
  std::vector<double> ap(n);

  double rsold = ParallelDotProduct(r, r);
  const double tolerance = 1e-8;

  // Индексы для итерации параллельного std::for_each
  std::vector<int> indices(n);
  std::iota(indices.begin(), indices.end(), 0);

  for (int iter = 0; iter < n * 2; ++iter) {
    // 1. Умножение матрицы на вектор (Matrix-Vector Multiply)
    std::for_each(std::execution::par, indices.begin(), indices.end(), [&](int i) {
      double sum = 0.0;
      size_t row_off = static_cast<size_t>(i) * n;
      for (int j = 0; j < n; ++j) {
        sum += a[row_off + j] * p[j];
      }
      ap[i] = sum;
    });

    double p_ap = ParallelDotProduct(p, ap);
    if (std::abs(p_ap) < 1e-15) {
      break;
    }

    const double alpha = rsold / p_ap;

    // 2. Обновление векторов x и r
    std::for_each(std::execution::par, indices.begin(), indices.end(), [&](int i) {
      x[i] += alpha * p[i];
      r[i] -= alpha * ap[i];
    });

    // КОРРЕКЦИЯ: Периодически пересчитываем истинную невязку r = b - A*x
    if ((iter + 1) % 10 == 0) {
      std::for_each(std::execution::par, indices.begin(), indices.end(), [&](int i) {
        double sum = 0.0;
        size_t row_off = static_cast<size_t>(i) * n;
        for (int j = 0; j < n; ++j) {
          sum += a[row_off + j] * x[j];
        }
        r[i] = b[i] - sum;
      });
    }

    const double rsnew = ParallelDotProduct(r, r);
    if (std::sqrt(rsnew) < tolerance) {
      break;
    }

    const double beta = rsnew / rsold;

    // 3. Обновление вектора p
    std::for_each(std::execution::par, indices.begin(), indices.end(), [&](int i) { p[i] = r[i] + beta * p[i]; });

    rsold = rsnew;
  }
  return true;
}

bool KruglovaAConjGradSleSTL::PostProcessingImpl() {
  return true;
}

}  // namespace kruglova_a_conjugate_gradient_sle
