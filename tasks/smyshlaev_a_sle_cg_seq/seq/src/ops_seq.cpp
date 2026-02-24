#include "smyshlaev_a_sle_cg_seq/seq/include/ops_seq.hpp"

#include <cmath>
#include <vector>

#include "smyshlaev_a_sle_cg_seq/common/include/common.hpp"

namespace smyshlaev_a_sle_cg_seq {

SmyshlaevASleCgTaskSEQ::SmyshlaevASleCgTaskSEQ(const InType &in) {
  SetTypeOfTask(GetStaticTypeOfTask());
  GetInput() = in;
}

bool SmyshlaevASleCgTaskSEQ::ValidationImpl() {
  const auto &A = GetInput().A;
  const auto &b = GetInput().b;
  if (A.empty() || b.empty()) {
    return false;
  }
  if (A.size() != b.size()) {
    return false;
  }
  if (A.size() != A[0].size()) {
    return false;
  }
  return true;
}

bool SmyshlaevASleCgTaskSEQ::PreProcessingImpl() {
  const auto &A = GetInput().A;
  int n = A.size();
  flat_A.resize(n * n);
  for (int i = 0; i < n; ++i) {
    for (int j = 0; j < n; ++j) {
      flat_A[i * n + j] = A[i][j];
    }
  }

  return true;
}

bool SmyshlaevASleCgTaskSEQ::RunImpl() {
  const auto &b = GetInput().b;
  int n = b.size();

  std::vector<double> r = b;
  std::vector<double> p = r;
  std::vector<double> Ap(n, 0.0);
  std::vector<double> result(n, 0.0);

  double rs_old = 0.0;
  for (int i = 0; i < n; ++i) {
    rs_old += r[i] * r[i];
  }

  const int max_iterations = n * 2;
  const double epsilon = 1e-9;

  // Если сразу всё хорошо — отдаём результат и выходим (не забываем GetOutput)
  if (std::sqrt(rs_old) < epsilon) {
    GetOutput() = result;
    return true;
  }

  for (int iter = 0; iter < max_iterations; ++iter) {
    // Ap = A * p
    for (int i = 0; i < n; ++i) {
      Ap[i] = 0.0;
      for (int j = 0; j < n; ++j) {
        Ap[i] += flat_A[i * n + j] * p[j];
      }
    }

    // Скалярное произведение p * Ap
    double pAp = 0.0;
    for (int i = 0; i < n; ++i) {
      pAp += p[i] * Ap[i];
    }

    // Защита от деления на ноль на случай плохой матрицы
    if (std::abs(pAp) < 1e-15) {
      break;
    }

    double alpha = rs_old / pAp;
    double rs_new = 0.0;

    // Обновляем вектора (без лишних аллокаций памяти)
    for (int i = 0; i < n; ++i) {
      result[i] += alpha * p[i];
      r[i] -= alpha * Ap[i];
      rs_new += r[i] * r[i];
    }

    if (std::sqrt(rs_new) < epsilon) {
      break;
    }

    double beta = rs_new / rs_old;
    for (int i = 0; i < n; ++i) {
      p[i] = r[i] + beta * p[i];
    }

    rs_old = rs_new;
  }

  GetOutput() = result;
  return true;
}

bool SmyshlaevASleCgTaskSEQ::PostProcessingImpl() {
  return true;
}

}  // namespace smyshlaev_a_sle_cg_seq