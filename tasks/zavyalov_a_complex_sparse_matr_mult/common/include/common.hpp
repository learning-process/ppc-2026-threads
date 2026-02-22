#pragma once

#include <string>
#include <tuple>

#include "task/include/task.hpp"

namespace zavyalov_a_compl_sparse_matr_mult {


    
using InType = int;
using OutType = int;
using TestType = std::tuple<int, std::string>;
using BaseTask = ppc::task::Task<InType, OutType>;

}  // namespace zavyalov_a_compl_sparse_matr_mult
