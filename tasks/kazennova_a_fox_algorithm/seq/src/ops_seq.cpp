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

  const int M = in.A.rows;  // число строк A
  const int N = in.B.cols;  // число столбцов B
  const int K = in.A.cols;  // общая размерность (столбцы A / строки B)

  const auto &A = in.A.data;
  const auto &B = in.B.data;
  auto &C = out.data;

  // Классическое умножение матриц: C[i][j] = sum_{k=0}^{K-1} A[i][k] * B[k][j]
  for (int i = 0; i < M; ++i) {
    for (int j = 0; j < N; ++j) {
      double sum = 0.0;
      for (int k = 0; k < K; ++k) {
        sum += A[static_cast<size_t>(i) * K + k] * B[static_cast<size_t>(k) * N + j];
      }
      C[static_cast<size_t>(i) * N + j] = sum;
    }
  }

  return true;
}

bool KazennovaATestTaskSEQ::PostProcessingImpl() {
  return !GetOutput().data.empty();
}

}  // namespace kazennova_a_fox_algorithm
