#include "kazennova_a_fox_algorithm/seq/include/ops_seq.hpp"

#include <cstddef>
#include <vector>

#include "kazennova_a_fox_algorithm/common/include/common.hpp"

namespace kazennova_a_fox_algorithm {

KazennovaATestTaskSEQ::KazennovaATestTaskSEQ(const InType &in) {
  SetTypeOfTask(GetStaticTypeOfTask());
  GetInput() = in;
}

bool KazennovaATestTaskSEQ::ValidationImpl() {
  const auto &in = GetInput();

  // Проверка на пустые данные
  if (in.A.data.empty() || in.B.data.empty()) {
    return false;
  }

  // Проверка размеров
  if (in.A.rows <= 0 || in.A.cols <= 0 || in.B.rows <= 0 || in.B.cols <= 0) {
    return false;
  }

  // Проверка согласованности матриц (A.cols == B.rows)
  if (in.A.cols != in.B.rows) {
    return false;
  }

  return true;
}

bool KazennovaATestTaskSEQ::PreProcessingImpl() {
  const auto &in = GetInput();

  // Выделяем память под результирующую матрицу C (размер A.rows × B.cols)
  GetOutput().rows = in.A.rows;
  GetOutput().cols = in.B.cols;
  GetOutput().data.assign(static_cast<size_t>(in.A.rows) * static_cast<size_t>(in.B.cols), 0.0);

  return true;
}

bool KazennovaATestTaskSEQ::RunImpl() {
  const auto &in = GetInput();
  auto &out = GetOutput();

  const int m = in.A.rows;
  const int n = in.B.cols;
  const int k = in.A.cols;

  const auto &a = in.A.data;
  const auto &b = in.B.data;
  auto &c = out.data;

  for (int i = 0; i < m; ++i) {
    for (int j = 0; j < n; ++j) {
      double sum = 0.0;
      for (int k_idx = 0; k_idx < k; ++k_idx) {
        const size_t a_idx = static_cast<size_t>(i) * k + k_idx;
        const size_t b_idx = static_cast<size_t>(k_idx) * n + j;
        sum += a[a_idx] * b[b_idx];
      }
      c[static_cast<size_t>(i) * n + j] = sum;
    }
  }

  return true;
}

bool KazennovaATestTaskSEQ::PostProcessingImpl() {
  return !GetOutput().data.empty();
}

}  // namespace kazennova_a_fox_algorithm
