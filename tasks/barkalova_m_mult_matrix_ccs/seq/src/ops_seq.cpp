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
/*
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

  // Базовая проверка (обязательно)
  if (A.rows <= 0 || A.cols <= 0 || B.rows <= 0 || B.cols <= 0 || A.cols != B.rows) {
    return false;
  }

  // Проверка структуры (критично для безопасности)
  if (A.col_ptrs.size() != static_cast<size_t>(A.cols) + 1 || 
      B.col_ptrs.size() != static_cast<size_t>(B.cols) + 1) {
    return false;
  }

  // Проверка целостности данных
  if (A.row_indices.size() != A.values.size() || 
      B.row_indices.size() != B.values.size()) {
    return false;
  }

  // Проверка базовых указателей
  if (A.col_ptrs.empty() || A.col_ptrs[0] != 0 || 
      B.col_ptrs.empty() || B.col_ptrs[0] != 0) {
    return false;
  }

  // Проверка индексов строк (чтобы избежать выхода за границы)
  for (int idx : A.row_indices) {
    if (idx < 0 || idx >= A.rows) return false;
  }
  for (int idx : B.row_indices) {
    if (idx < 0 || idx >= B.rows) return false;
  }

  return true;
}

bool BarkalovaMMultMatrixCcsSEQ::PreProcessingImpl() {
  return true;
}
//это было у меня написано
/*
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
      if (pos >= dst.nnz) {
        continue;  // Проверка границ
      }

      dst.values[pos] = val;
      dst.row_indices[pos] = col;
      current_pos[row]++;
    }
  }
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
  // (это будут столбцы транспонированной матрицы)
  std::vector<int> row_count(dst.cols, 0);  // dst.cols = src.rows
  for (int i = 0; i < src.nnz; ++i) {
    int row = src.row_indices[i];
    row_count[row]++;
  }

  // Построение col_ptrs для транспонированной матрицы
  dst.col_ptrs.resize(dst.cols + 1);
  dst.col_ptrs[0] = 0;
  for (int i = 0; i < dst.cols; ++i) {
    dst.col_ptrs[i + 1] = dst.col_ptrs[i] + row_count[i];
  }

  // Временный массив для отслеживания текущих позиций
  std::vector<int> current_pos(dst.cols, 0);
  dst.values.resize(src.nnz);
  dst.row_indices.resize(src.nnz);

  // Заполнение values и row_indices для транспонированной матрицы
  for (int col = 0; col < src.cols; ++col) {
    for (int ip = src.col_ptrs[col]; ip < src.col_ptrs[col + 1]; ++ip) {
      int row = src.row_indices[ip];
      Complex val = src.values[ip];

      int pos = dst.col_ptrs[row] + current_pos[row];
      dst.values[pos] = val;
      dst.row_indices[pos] = col;  // col становится индексом строки в транспонированной
      current_pos[row]++;
    }
  }
}
//изначально было написано вот так, но вроде есть ошибка, у власовой через маркеры сделано
/*
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






void BarkalovaMMultMatrixCcsSEQ::MultiplyMatrices(const CCSMatrix &a, const CCSMatrix &b, CCSMatrix &c) {
  // Транспонируем A для эффективного доступа к строкам
  CCSMatrix at;
  TransposeMatrix(a, at);

  // Инициализация результата
  c.rows = a.rows;
  c.cols = b.cols;
  c.col_ptrs.clear();
  c.col_ptrs.push_back(0);
  c.values.clear();
  c.row_indices.clear();

  // Для каждого столбца j матрицы B
  for (int j = 0; j < b.cols; ++j) {
    // Временный вектор для текущего столбца результата
    std::vector<Complex> column(c.rows, Complex(0.0, 0.0));

    // Для каждого ненулевого элемента в столбце j матрицы B
    for (int kp = b.col_ptrs[j]; kp < b.col_ptrs[j + 1]; ++kp) {
      int k = b.row_indices[kp];  // индекс строки в B
      Complex b_val = b.values[kp];

      // В транспонированной матрице at столбец k содержит элементы строки k исходной A
      for (int ip = at.col_ptrs[k]; ip < at.col_ptrs[k + 1]; ++ip) {
        int i = at.row_indices[ip];  // индекс строки в A
        Complex a_val = at.values[ip];
        
        // Накапливаем результат для позиции (i, j)
        column[i] += a_val * b_val;
      }
    }

    // Сохраняем ненулевые элементы текущего столбца
    for (int i = 0; i < c.rows; ++i) {
      if (std::norm(column[i]) > kEpsilon) {
        c.values.push_back(column[i]);
        c.row_indices.push_back(i);
      }
    }

    c.col_ptrs.push_back(static_cast<int>(c.values.size()));
  }

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
*/




//как у власовой
#include "barkalova_m_mult_matrix_ccs/seq/include/ops_seq.hpp"

#include <cmath>
#include <cstddef>
#include <exception>
#include <vector>
#include <complex>

#include "barkalova_m_mult_matrix_ccs/common/include/common.hpp"

namespace barkalova_m_mult_matrix_ccs {

namespace {
constexpr double kEpsilon = 1e-10;
}  // namespace

BarkalovaMMultMatrixCcsSEQ::BarkalovaMMultMatrixCcsSEQ(const InType &in) {
  SetTypeOfTask(GetStaticTypeOfTask());
  GetInput() = in;
  GetOutput() = CCSMatrix{};
}

bool BarkalovaMMultMatrixCcsSEQ::ValidationImpl() {
  const auto &[A, B] = GetInput();

  // Базовая проверка размерностей
  if (A.rows <= 0 || A.cols <= 0 || B.rows <= 0 || B.cols <= 0 || A.cols != B.rows) {
    return false;
  }

  // Проверка структуры CCS
  if (A.col_ptrs.size() != static_cast<size_t>(A.cols) + 1 || 
      B.col_ptrs.size() != static_cast<size_t>(B.cols) + 1) {
    return false;
  }

  // Проверка целостности данных
  if (A.row_indices.size() != A.values.size() || 
      B.row_indices.size() != B.values.size()) {
    return false;
  }

  // Проверка, что col_ptrs не пустые и начинаются с 0
  if (A.col_ptrs.empty() || A.col_ptrs[0] != 0 || 
      B.col_ptrs.empty() || B.col_ptrs[0] != 0) {
    return false;
  }

  // Проверка, что col_ptrs монотонно возрастают
  for (size_t i = 0; i < A.col_ptrs.size() - 1; ++i) {
    if (A.col_ptrs[i] > A.col_ptrs[i + 1]) return false;
  }
  for (size_t i = 0; i < B.col_ptrs.size() - 1; ++i) {
    if (B.col_ptrs[i] > B.col_ptrs[i + 1]) return false;
  }

  // Проверка индексов строк (только если есть ненулевые элементы)
  for (int idx : A.row_indices) {
    if (idx < 0 || idx >= A.rows) return false;
  }
  for (int idx : B.row_indices) {
    if (idx < 0 || idx >= B.rows) return false;
  }

  // Проверка, что количество ненулевых элементов соответствует col_ptrs
  if (A.nnz != static_cast<int>(A.values.size()) || 
      B.nnz != static_cast<int>(B.values.size())) {
    return false;
  }
  
  if (A.nnz != A.col_ptrs.back() || B.nnz != B.col_ptrs.back()) {
    return false;
  }

  return true;
}

bool BarkalovaMMultMatrixCcsSEQ::PreProcessingImpl() {
  return true;
}

void BarkalovaMMultMatrixCcsSEQ::TransposeMatrix(const CCSMatrix &a, CCSMatrix &at) {
  at.rows = a.cols;
  at.cols = a.rows;
  at.nnz = a.nnz;

  if (a.nnz == 0) {
    at.values.clear();
    at.row_indices.clear();
    at.col_ptrs.assign(at.cols + 1, 0);
    return;
  }

  std::vector<int> row_count(at.cols, 0);
  for (int i = 0; i < a.nnz; i++) {
    row_count[a.row_indices[i]]++;
  }

  at.col_ptrs.resize(at.cols + 1);
  at.col_ptrs[0] = 0;
  for (int i = 0; i < at.cols; i++) {
    at.col_ptrs[i + 1] = at.col_ptrs[i] + row_count[i];
  }

  at.values.resize(a.nnz);
  at.row_indices.resize(a.nnz);

  std::vector<int> current_pos(at.cols, 0); 
  for (int col = 0; col < a.cols; col++) {
    for (int i = a.col_ptrs[col]; i < a.col_ptrs[col + 1]; i++) {
      int row = a.row_indices[i];
      Complex val = a.values[i];

      int pos = at.col_ptrs[row] + current_pos[row];
      at.values[pos] = val;
      at.row_indices[pos] = col;
      current_pos[row]++;
    }
  }
}

namespace {

void ProcessColumnSEQ(const CCSMatrix &at, const CCSMatrix &b, int col_index, 
                      std::vector<Complex> &temp_row, std::vector<int> &row_marker,
                      std::vector<Complex> &res_val, std::vector<int> &res_row_ind) {
  for (int k = b.col_ptrs[col_index]; k < b.col_ptrs[col_index + 1]; k++) {
    int row_b = b.row_indices[k];
    Complex val_b = b.values[k];

    for (int idx = at.col_ptrs[row_b]; idx < at.col_ptrs[row_b + 1]; idx++) {
      int row_a = at.row_indices[idx];
      Complex val_a = at.values[idx];

      if (row_marker[row_a] != col_index) {
        row_marker[row_a] = col_index;
        temp_row[row_a] = val_a * val_b;
      } else {
        temp_row[row_a] += val_a * val_b;
      }
    }
  }

  for (size_t i = 0; i < temp_row.size(); i++) {
    if (row_marker[i] == col_index && std::norm(temp_row[i]) > kEpsilon) {
      res_val.push_back(temp_row[i]);
      res_row_ind.push_back(static_cast<int>(i));
    }
  }
}

}  // namespace

void BarkalovaMMultMatrixCcsSEQ::MultiplyMatrices(const CCSMatrix &a, const CCSMatrix &b, CCSMatrix &c) {
  CCSMatrix at;
  TransposeMatrix(a, at);

  c.rows = a.rows;
  c.cols = b.cols;
  c.col_ptrs.push_back(0);

  std::vector<Complex> temp_row(c.rows, Complex(0.0, 0.0));
  std::vector<int> row_marker(c.rows, -1);

  for (int j = 0; j < b.cols; j++) {
    ProcessColumnSEQ(at, b, j, temp_row, row_marker, c.values, c.row_indices);
    c.col_ptrs.push_back(static_cast<int>(c.values.size()));
  }

  c.nnz = static_cast<int>(c.values.size());
}

bool BarkalovaMMultMatrixCcsSEQ::RunImpl() {
  const auto &[a, b] = GetInput();

  try {
    CCSMatrix c;
    MultiplyMatrices(a, b, c);
    GetOutput() = c;
    return true;
  } catch (const std::exception &) {
    return false;
  }
}

bool BarkalovaMMultMatrixCcsSEQ::PostProcessingImpl() {
  const auto &c = GetOutput();
  return c.rows > 0 && c.cols > 0 && c.col_ptrs.size() == static_cast<size_t>(c.cols) + 1;
}

}  // namespace barkalova_m_mult_matrix_ccs