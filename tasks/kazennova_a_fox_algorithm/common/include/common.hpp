#pragma once

#include <cstddef>
#include <string>
#include <tuple>
#include <vector>

#include "task/include/task.hpp"

namespace kazennova_a_fox_algorithm {

struct Matrix {
  std::vector<double> data;  // элементы по строкам
  int rows;
  int cols;
};

struct MatricesPair {
  Matrix A;
  Matrix B;
};

using InType = MatricesPair;
using OutType = Matrix;
using TestType = std::tuple<int, std::string>;
using BaseTask = ppc::task::Task<InType, OutType>;

}  // namespace kazennova_a_fox_algorithm
