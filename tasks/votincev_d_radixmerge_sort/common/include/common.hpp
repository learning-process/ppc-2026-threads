#pragma once

#include <string>
#include <vector>

#include "task/include/task.hpp"

namespace votincev_d_radixmerge_sort {

using InType = std::vector<int>;
using OutType = std::vector<int>;
using TestType = std::string;
using BaseTask = ppc::task::Task<InType, OutType>;

}  // namespace votincev_d_radixmerge_sort
