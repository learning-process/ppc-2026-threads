#pragma once

#include <cstddef>
#include <string>
#include <tuple>
#include <vector>

#include "task/include/task.hpp"

namespace kulik_a_mat_mul_double_ccs {

struct CCS {
  size_t n = 0;   // кол-во строк
  size_t m = 0;   // кол-во столбцов
  size_t nz = 0;  // кол-во ненулевых элементов
  std::vector<size_t> col_ind;
  std::vector<size_t> row;
  std::vector<double> value;

  CCS() = default;
  CCS(const CCS &) = default;
  CCS &operator=(const CCS &) = default;
};

using InType = std::tuple<CCS, CCS>;
using OutType = CCS;
using TestType = std::tuple<std::string, std::string>;
using BaseTask = ppc::task::Task<InType, OutType>;

}  // namespace kulik_a_mat_mul_double_ccs
