/*#include "barkalova_m_mult_matrix_ccs/seq/include/ops_seq.hpp"

#include <numeric>
#include <vector>

#include "barkalova_m_mult_matrix_ccs/common/include/common.hpp"
#include "util/include/util.hpp"

namespace barkalova_m_mult_matrix_ccs {

BarkalobaMMultMatrixCcsSEQ::BarkalobaMMultMatrixCcsSEQ(const InType &in) {
  SetTypeOfTask(GetStaticTypeOfTask());
  GetInput() = in;
  GetOutput() = 0;
}

bool BarkalobaMMultMatrixCcsSEQ::ValidationImpl() {
  return (GetInput() > 0) && (GetOutput() == 0);
}

bool BarkalobaMMultMatrixCcsSEQ::PreProcessingImpl() {
  GetOutput() = 2 * GetInput();
  return GetOutput() > 0;
}

bool BarkalobaMMultMatrixCcsSEQ::RunImpl() {
  if (GetInput() == 0) {
    return false;
  }

  for (InType i = 0; i < GetInput(); i++) {
    for (InType j = 0; j < GetInput(); j++) {
      for (InType k = 0; k < GetInput(); k++) {
        std::vector<InType> tmp(i + j + k, 1);
        GetOutput() += std::accumulate(tmp.begin(), tmp.end(), 0);
        GetOutput() -= i + j + k;
      }
    }
  }

  const int num_threads = ppc::util::GetNumThreads();
  GetOutput() *= num_threads;

  int counter = 0;
  for (int i = 0; i < num_threads; i++) {
    counter++;
  }

  if (counter != 0) {
    GetOutput() /= counter;
  }
  return GetOutput() > 0;
}

bool BarkalobaMMultMatrixCcsSEQ::PostProcessingImpl() {
  GetOutput() -= GetInput();
  return GetOutput() > 0;
}

}  // namespace barkalova_m_mult_matrix_ccs
*/


#include "barkalova_m_mult_matrix_ccs/seq/include/ops_seq.hpp"

#include <algorithm>
#include <cmath>
#include <complex>
#include <exception>
#include <vector>

namespace barkalova_m_mult_matrix_ccs {

BarkalovaMMultMatrixCcsSEQ::BarkalovaMMultMatrixCcsSEQ(const InType &in) {
  SetTypeOfTask(GetStaticTypeOfTask());
  GetInput() = in;
  GetOutput() = CCSMatrix();
}

bool BarkalovaMMultMatrixCcsSEQ::ValidationImpl() {
  const auto &[A, B] = GetInput();

  // Проверка размеров матриц
  if (A.rows <= 0 || A.cols <= 0 || B.rows <= 0 || B.cols <= 0) {
    return false;
  }

  // Проверка согласованности размерностей для умножения
  if (A.cols != B.rows) {
    return false;
  }

  // Проверка корректности структуры CCS
  if (A.col_ptrs.size() != static_cast<size_t>(A.cols) + 1 || B.col_ptrs.size() != static_cast<size_t>(B.cols) + 1) {
    return false;
  }

  if (A.row_indices.size() != A.values.size() || B.row_indices.size() != B.values.size()) {
    return false;
  }

  // Проверка, что первый указатель столбцов равен 0
  if ((!A.col_ptrs.empty() && A.col_ptrs[0] != 0) || (!B.col_ptrs.empty() && B.col_ptrs[0] != 0)) {
    return false;
  }

  // Проверка корректности указателей столбцов (неубывающая последовательность)
  for (size_t j = 0; j < A.col_ptrs.size() - 1; ++j) {
    if (A.col_ptrs[j] > A.col_ptrs[j + 1]) {
      return false;
    }
  }

  for (size_t j = 0; j < B.col_ptrs.size() - 1; ++j) {
    if (B.col_ptrs[j] > B.col_ptrs[j + 1]) {
      return false;
    }
  }

  // Проверка индексов строк в допустимых пределах
  for (size_t i = 0; i < A.row_indices.size(); ++i) {
    if (A.row_indices[i] < 0 || A.row_indices[i] >= A.rows) {
      return false;
    }
  }

  for (size_t i = 0; i < B.row_indices.size(); ++i) {
    if (B.row_indices[i] < 0 || B.row_indices[i] >= B.rows) {
      return false;
    }
  }

  // Проверка упорядоченности индексов строк в каждом столбце
  for (int j = 0; j < A.cols; ++j) {
    for (int p = A.col_ptrs[j]; p < A.col_ptrs[j + 1] - 1; ++p) {
      if (A.row_indices[p] >= A.row_indices[p + 1]) {
        return false;
      }
    }
  }

  for (int j = 0; j < B.cols; ++j) {
    for (int p = B.col_ptrs[j]; p < B.col_ptrs[j + 1] - 1; ++p) {
      if (B.row_indices[p] >= B.row_indices[p + 1]) {
        return false;
      }
    }
  }

  return true;
}

bool BarkalovaMMultMatrixCcsSEQ::PreProcessingImpl() {
  // Очищаем результат перед вычислениями
  GetOutput() = CCSMatrix();  // так ли и обязаткльно ли это надо писать
  return true;
}

void BarkalovaMMultMatrixCcsSEQ::TransposeMatrix(const CCSMatrix &src, CCSMatrix &dst) {
  // Устанавливаем размеры транспонированной матрицы
  dst.rows = src.cols;
  dst.cols = src.rows;
  dst.nnz = src.nnz;

  if (src.nnz == 0) {
    dst.col_ptrs.assign(dst.cols + 1, 0);
    dst.values.clear();
    dst.row_indices.clear();
    return;
  }

  // Подсчет количества элементов в каждой строке исходной матрицы
  std::vector<int> row_count(src.rows, 0);
  for (int i = 0; i < src.nnz; ++i) {
    int row = src.row_indices[i];
    if (row >= 0 && row < src.rows) {
      row_count[row]++;
    }
  }

  // Построение col_ptrs для транспонированной матрицы
  dst.col_ptrs.resize(dst.cols + 1);
  dst.col_ptrs[0] = 0;
  for (int i = 0; i < dst.cols; ++i) {
    dst.col_ptrs[i + 1] = dst.col_ptrs[i] + row_count[i];
  }

  // Временный массив для отслеживания текущих позиций
  std::vector<int> current_pos(dst.cols, 0);

  // Заполнение values и row_indices для транспонированной матрицы
  dst.values.resize(src.nnz);
  dst.row_indices.resize(src.nnz);

  for (int col = 0; col < src.cols; ++col) {
    for (int ip = src.col_ptrs[col]; ip < src.col_ptrs[col + 1]; ++ip) {
      int row = src.row_indices[ip];
      if (row < 0 || row >= dst.cols) {
        continue;  // Проверка границ
      }

      Complex val = src.values[ip];

      int pos = dst.col_ptrs[row] + current_pos[row];
      if (pos >= src.nnz) {
        continue;  // Проверка границ
      }

      dst.values[pos] = val;
      dst.row_indices[pos] = col;
      current_pos[row]++;
    }
  }
}

void BarkalovaMMultMatrixCcsSEQ::MultiplyMatrices(const CCSMatrix &a, const CCSMatrix &b, CCSMatrix &c) {
  // Транспонируем первую матрицу
  CCSMatrix at;
  TransposeMatrix(a, at);

  // Инициализация результирующей матрицы
  c.rows = a.rows;
  c.cols = b.cols;
  c.col_ptrs.assign(c.cols + 1, 0);
  c.values.clear();
  c.row_indices.clear();

  // Для каждого столбца j матрицы B
  for (int j = 0; j < b.cols; ++j) {
    // Временный вектор для хранения результатов текущего столбца
    std::vector<Complex> column(c.rows, Complex(0.0, 0.0));

    // Для каждого ненулевого элемента в столбце j матрицы B
    for (int kp = b.col_ptrs[j]; kp < b.col_ptrs[j + 1]; ++kp) {
      int k = b.row_indices[kp];  // индекс строки в B
      if (k < 0 || k >= a.cols) {
        continue;  // Проверка границ
      }

      Complex b_val = b.values[kp];

      // В транспонированной матрице at столбец k соответствует строке k исходной A
      for (int ip = at.col_ptrs[k]; ip < at.col_ptrs[k + 1]; ++ip) {
        int i = at.row_indices[ip];  // индекс строки в A
        if (i < 0 || i >= c.rows) {
          continue;  // Проверка границ
        }

        Complex a_val = at.values[ip];
        column[i] += a_val * b_val;
      }
    }

    // Запоминаем начало столбца j
    c.col_ptrs[j] = static_cast<int>(c.values.size());

    // Сохраняем ненулевые элементы столбца
    for (int i = 0; i < c.rows; ++i) {
      if (std::norm(column[i]) > kEpsilon) {
        c.values.push_back(column[i]);
        c.row_indices.push_back(i);
      }
    }
  }

  // Устанавливаем последний указатель
  c.col_ptrs[c.cols] = static_cast<int>(c.values.size());
  c.nnz = static_cast<int>(c.values.size());
}

bool BarkalovaMMultMatrixCcsSEQ::RunImpl() {
  const auto &[A, B] = GetInput();

  try {
    CCSMatrix result;
    MultiplyMatrices(A, B, result);
    GetOutput() = result;
    return true;
  } catch (const std::exception &) {  // Убрали имя переменной
    return false;
  }
}

bool BarkalovaMMultMatrixCcsSEQ::PostProcessingImpl() {
  // Простая проверка результата
  const auto &C = GetOutput();
  const auto &[A, B] = GetInput();

  if (C.rows != A.rows || C.cols != B.cols) {
    return false;
  }

  return true;
}

}  // namespace barkalova_m_mult_matrix_ccs