#pragma once

#include <complex>
#include <cstddef>
#include <string>
#include <tuple>
#include <utility>
#include <vector>

#include "task/include/task.hpp"

namespace posternak_a_crs_mul_complex_matrix {

struct CRSMatrix {
  int rows = 0;
  int cols = 0;
  std::vector<std::complex<double>> values;  // вектор ненулевых значений
  std::vector<int> index_col;                // индексация столбцов для каждого значения
  std::vector<int> index_row;  // количество ненулевых элементов в строках (размер N+1 для запирающего 0 в начале)

  [[nodiscard]] bool IsValid() const {
    if (rows < 0 || cols < 0) {
      return false;
    }
    if (index_row.size() != static_cast<size_t>(rows) + 1) {
      return false;
    }
    if (values.size() != index_col.size()) {
      return false;
    }
    if (!values.empty() && std::cmp_not_equal(index_row.back(), static_cast<int>(values.size()))) {
      return false;
    }
    return true;
  }
};

using InType = std::pair<CRSMatrix, CRSMatrix>;
using OutType = CRSMatrix;
using TestType = std::tuple<int, std::string>;
using BaseTask = ppc::task::Task<InType, OutType>;

CRSMatrix MakeCRS(int rows, int cols, const std::vector<std::vector<std::complex<double>>> &default_matrix);
bool MatricesEqual(const CRSMatrix &a, const CRSMatrix &b);

CRSMatrix MakeBandedCRS(int size, int bandwidth);
std::complex<double> ComputeExpectedValue(const CRSMatrix &a, const CRSMatrix &b, int row, int col);
bool CheckKeyElements(const CRSMatrix &result, const CRSMatrix &a, const CRSMatrix &b);
bool CheckSingleElement(const CRSMatrix &result, const CRSMatrix &a, const CRSMatrix &b, int row, int col);

}  // namespace posternak_a_crs_mul_complex_matrix
