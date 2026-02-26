#pragma once

#include <string>
#include <tuple>
#include <complex>

#include "task/include/task.hpp"

namespace liulin_y_complex_ccs {
struct CCSMatrix{
    std::vector<std::complex<double>> values;
    std::vector<int> col_index;
    std::vector<int> row_index;
    int count_rows = 0;
    int count_cols = 0;
};
using InType = std::pair<CCSMatrix, CCSMatrix>;
using OutType = CCSMatrix;
using TestType = std::tuple<int, std::string>;
using BaseTask = ppc::task::Task<InType, OutType>;

}  // namespace liulin_y_complex_ccs
