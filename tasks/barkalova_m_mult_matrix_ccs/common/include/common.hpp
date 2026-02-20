/*#pragma once

#include <string>
#include <tuple>

#include "task/include/task.hpp"

namespace barkalova_m_mult_matrix_ccs {

using InType = int;
using OutType = int;
using TestType = std::tuple<int, std::string>;
using BaseTask = ppc::task::Task<InType, OutType>;

}  // namespace barkalova_m_mult_matrix_ccs
*/

#pragma once

#include <complex>
#include <cstddef>
#include <tuple>
#include <utility>
#include <vector>

#include "task/include/task.hpp"

namespace barkalova_m_mult_matrix_ccs {

using Complex = std::complex<double>;

// Структура для разреженной матрицы в формате CCS
struct CCSMatrix {
  std::vector<Complex> values;   // ненулевые значения
  std::vector<int> row_indices;  // индексы строк для каждого значения
  std::vector<int> col_ptrs;     // указатели на начало каждого столбца
  int rows = 0;                  // количество строк
  int cols = 0;                  // количество столбцов
  int nnz = 0;                   // количество ненулевых элементов

  CCSMatrix() = default;

  CCSMatrix(int r, int c) : rows(r), cols(c) {
    col_ptrs.resize(cols + 1, 0);
  }

  // Очистка матрицы
  void Clear() {
    values.clear();
    row_indices.clear();
    col_ptrs.assign(cols + 1, 0);
    nnz = 0;
  }

  // Проверка на пустоту
  bool Empty() const {
    return values.empty();
  }

  // Количество ненулевых элементов - ИМЕННО ТАК, КАК В ТЕСТАХ!
  size_t nonZeros() const {
    return values.size();
  }
};

using InType = std::pair<CCSMatrix, CCSMatrix>;
using OutType = CCSMatrix;
using TestType = std::tuple<int, int, int>;  // rows, inner_dim, cols
using BaseTask = ppc::task::Task<InType, OutType>;

}  // namespace barkalova_m_mult_matrix_ccs
