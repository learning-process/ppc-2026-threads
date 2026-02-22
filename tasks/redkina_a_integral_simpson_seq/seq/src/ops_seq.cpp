#include "redkina_a_integral_simpson_seq/seq/include/ops_seq.hpp"

#include <cmath>
#include <cstddef>
#include <functional>
#include <vector>

#include "redkina_a_integral_simpson_seq/common/include/common.hpp"

namespace redkina_a_integral_simpson_seq {

RedkinaAIntegralSimpsonSEQ::RedkinaAIntegralSimpsonSEQ(const InType &in) {
  SetTypeOfTask(GetStaticTypeOfTask());
  GetInput() = in;
}

bool RedkinaAIntegralSimpsonSEQ::ValidationImpl() {
  const auto &in = GetInput();
  size_t dim = in.a.size();

  if (dim == 0 || in.b.size() != dim || in.n.size() != dim) {
    return false;
  }

  for (size_t i = 0; i < dim; ++i) {
    if (in.a[i] >= in.b[i]) {
      return false;
    }
    if (in.n[i] <= 0 || in.n[i] % 2 != 0) {
      return false;
    }
  }

  if (!in.func) {
    return false;
  }

  return true;
}

bool RedkinaAIntegralSimpsonSEQ::PreProcessingImpl() {
  const auto &in = GetInput();
  func_ = in.func;
  a_ = in.a;
  b_ = in.b;
  n_ = in.n;
  result_ = 0.0;
  return true;
}

bool RedkinaAIntegralSimpsonSEQ::RunImpl() {
  size_t dim = a_.size();

  // Шаги по каждому измерению
  std::vector<double> h(dim);
  for (size_t i = 0; i < dim; ++i) {
    h[i] = (b_[i] - a_[i]) / static_cast<double>(n_[i]);
  }

  // Произведение шагов (h1 * h2 * ... * hd)
  double h_prod = 1.0;
  for (size_t i = 0; i < dim; ++i) {
    h_prod *= h[i];
  }

  std::vector<double> point(dim);
  double sum = 0.0;

  // Массив текущих индексов по каждому измерению
  std::vector<int> indices(dim, 0);

  // Итеративный перебор всех узлов сетки (многомерный счётчик)
  while (true) {
    // Вычисление координат точки и произведения весов Симпсона
    double w_prod = 1.0;
    for (size_t d = 0; d < dim; ++d) {
      int idx = indices[d];
      point[d] = a_[d] + static_cast<double>(idx) * h[d];

      // Вес Симпсона для данного индекса
      int w = 0;
      if (idx == 0 || idx == n_[d]) {
        w = 1;
      } else if (idx % 2 == 1) {
        w = 4;
      } else {
        w = 2;
      }
      w_prod *= static_cast<double>(w);
    }

    sum += w_prod * func_(point);

    // Переход к следующей комбинации индексов (как в счётчике с переменным основанием)
    int d = static_cast<int>(dim) - 1;
    while (d >= 0 && indices[d] == n_[d]) {
      indices[d] = 0;
      --d;
    }
    if (d < 0) {
      break;  // все комбинации перебраны
    }
    ++indices[d];
  }

  // Знаменатель: 3^dim
  double denominator = 1.0;
  for (size_t i = 0; i < dim; ++i) {
    denominator *= 3.0;
  }

  result_ = (h_prod / denominator) * sum;
  return true;
}

bool RedkinaAIntegralSimpsonSEQ::PostProcessingImpl() {
  GetOutput() = result_;
  return true;
}

}  // namespace redkina_a_integral_simpson_seq
