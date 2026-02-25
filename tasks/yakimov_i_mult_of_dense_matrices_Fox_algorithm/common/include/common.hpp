#pragma once

#include <string>
#include <tuple>
#include <vector>

#include "task/include/task.hpp"

namespace yakimov_i_mult_of_dense_matrices_Fox_algorithm {

using InType = int;
using OutType = double;
using TestType = std::tuple<int, std::string>;
using BaseTask = ppc::task::Task<InType, OutType>;

struct DenseMatrix {
    std::vector<double> data;
    int rows = 0;
    int cols = 0;
    
    DenseMatrix() = default;
    
    double& operator()(int i, int j) {
        return data[static_cast<size_t>(i * cols + j)];
    }
    
    const double& operator()(int i, int j) const {
        return data[static_cast<size_t>(i * cols + j)];
    }
};

}  // namespace yakimov_i_mult_of_dense_matrices_Fox_algorithm