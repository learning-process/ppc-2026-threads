#pragma once

#include <string>
#include <tuple>
#include <vector>

#include "task/include/task.hpp"

namespace marin_l_mark_components_seq {

using Image = std::vector<std::vector<int>>;
using Labels = std::vector<std::vector<int>>;

struct InData {
  Image binary;
};

struct OutData {
  Labels labels;
};

using InType = InData;
using OutType = OutData;
using TestType = std::tuple<int, int, std::string>;
using BaseTask = ppc::task::Task<InType, OutType>;

}  // namespace marin_l_mark_components_seq
