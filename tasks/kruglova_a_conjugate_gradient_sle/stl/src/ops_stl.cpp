#include "kruglova_a_conjugate_gradient_sle/stl/include/ops_stl.hpp"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <execution>
#include <numeric>
#include <vector>

#include "kruglova_a_conjugate_gradient_sle/common/include/common.hpp"

namespace kruglova_a_conjugate_gradient_sle {

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

  if (n <= 0) {
    return true;
  }

  std::vector<double> r = b;
  std::vector<double> p = r;
  std::vector<double> ap(n, 0.0);

  // Вектор индексов для организации параллельного обхода через std::for_each
  std::vector<int> indices(n);
  std::iota(indices.begin(), indices.end(), 0);

  // Вычисление стартового значения rsold через параллельный скалярный редуктор
  double rsold = std::transform_reduce(std::execution::par, r.begin(), r.end(), r.begin(), 0.0);

  const double tolerance = 1e-8;
  const int max_iter = n * 2;

  for (int iter = 0; iter < max_iter; ++iter) {
    // 1. Матрично-векторное умножение: ap = A * p
    std::for_each(std::execution::par, indices.begin(), indices.end(), [&](int i) {
      double sum = 0.0;
      size_t row_offset = static_cast<size_t>(i) * n;
      for (int j = 0; j < n; ++j) {
        sum += a[row_offset + j] * p[j];
      }
      ap[i] = sum;
    });

    // 2. Вычисление скалярного произведения: p_ap = p * ap
    double p_ap = std::transform_reduce(std::execution::par, p.begin(), p.end(), ap.begin(), 0.0);

    if (std::abs(p_ap) < 1e-15) {
      break;
    }

    const double alpha = rsold / p_ap;

    // 3. Параллельное обновление векторов x и r
    std::for_each(std::execution::par, indices.begin(), indices.end(), [&](int i) {
      x[i] += alpha * p[i];
      r[i] -= alpha * ap[i];
    });

    // 4. Расчет новой невязки rsnew = r * r
    double rsnew = std::transform_reduce(std::execution::par, r.begin(), r.end(), r.begin(), 0.0);

    if (std::sqrt(rsnew) < tolerance) {
      break;
    }

    const double beta = rsnew / rsold;

    // 5. Параллельное обновление вектора направлений p
    std::for_each(std::execution::par, indices.begin(), indices.end(), [&](int i) { p[i] = r[i] + beta * p[i]; });

    rsold = rsnew;
  }

  return true;
}

bool KruglovaAConjGradSleSTL::PostProcessingImpl() {
  return true;
}

}  // namespace kruglova_a_conjugate_gradient_sle
