#include "makoveeva_matmul_double_omp/omp/include/ops_omp.hpp"

#include <omp.h>

#include <cmath>
#include <cstddef>
#include <vector>

#include "makoveeva_matmul_double_omp/common/include/common.hpp"

namespace makoveeva_matmul_double_omp {

// Убираем : n_(0) - инициализация будет в классе
MatmulDoubleOMPTask::MatmulDoubleOMPTask(const InType &in) {
  SetTypeOfTask(GetStaticTypeOfTask());
  GetInput() = in;
  GetOutput() = std::vector<double>();
}

bool MatmulDoubleOMPTask::ValidationImpl() {
  const auto &input = GetInput();
  const size_t n = std::get<0>(input);
  const auto &a = std::get<1>(input);
  const auto &b = std::get<2>(input);

  return n > 0 && a.size() == n * n && b.size() == n * n;
}

bool MatmulDoubleOMPTask::PreProcessingImpl() {
  const auto &input = GetInput();
  n_ = std::get<0>(input);
  A_ = std::get<1>(input);
  B_ = std::get<2>(input);
  C_.assign(n_ * n_, 0.0);

  return true;
}

bool MatmulDoubleOMPTask::RunImpl() {
  if (n_ <= 0) {
    return false;
  }

  const size_t n = n_;
  const auto &a = A_;
  const auto &b = B_;
  auto &c = C_;

#pragma omp parallel for collapse(2) default(none) shared(a, b, c, n)
  for (size_t i = 0; i < n; ++i) {
    for (size_t j = 0; j < n; ++j) {
      double sum = 0.0;
      for (size_t k = 0; k < n; ++k) {
        sum += a[(i * n) + k] * b[(k * n) + j];
      }
      c[(i * n) + j] = sum;
    }
  }

  GetOutput() = C_;
  return true;
}

bool MatmulDoubleOMPTask::PostProcessingImpl() {
  return true;
}

}  // namespace makoveeva_matmul_double_omp
